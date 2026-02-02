from modeling.simulation.study_params import StudyParams
from modeling.models import model_zoo
from modeling.models.model_iterface import SimulationParams
from modeling.models.model_zoo import ModelNames
from modeling.simulation.simulations_db import SimulationsDB
from datetime import datetime as dt
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from typing import Any, Tuple, List
from tqdm import tqdm
from itertools import product

DATABASE_PATH = "data/comb_simulations"


def _generate_datasets_for_params(
    args: Tuple[SimulationParams, int, ModelNames, str],
) -> Tuple[SimulationParams, List[Any]]:
    """
    Worker function: Generate datasets for a single parameter combination.

    Each process gets its own model instance.
    Returns sim_params and list of generated datasets.
    """
    sim_params, num_datasets, sim_name_for_data, _ = args

    # Each worker creates its own model instance
    model = model_zoo.get_model_by_name(sim_name_for_data)

    datasets = []
    for _ in range(num_datasets):
        data = model.simulate_data(sim_params)
        datasets.append(data)

    return sim_params, datasets


def make_datasets(
    db: SimulationsDB,
    params: StudyParams,
    name: str,
    num_datasets: int,
    num_processes: int | None = None,
) -> int:
    """
    Create study and generate datasets for all parameter combinations in parallel.

    Args:
        db: Database connection
        params: Study parameters
        name: Study name
        num_datasets: Number of datasets to generate per parameter combination
        num_processes: Number of parallel processes (defaults to CPU count)
    """
    if num_processes is None:
        num_processes = cpu_count()

    study_id = db.insert_study(name, dt.now(), params)

    param_combinations = product(
        params.beta0,
        params.beta1,
        params.p11,
        params.p_aru11,
        params.p_aru01,
        params.mu,
        params.sigma,
        params.n_sites,
        params.n_surveys_pc,
        params.n_surveys_aru,
        params.aru_scores_independent_data,
        params.n_surveys_scores,
    )

    # Build list of valid parameter combinations
    sim_args = []
    for (
        beta0,
        beta1,
        p11,
        p_aru11,
        p_aru01,
        mu,
        sigma,
        n_sites,
        n_surveys_pc,
        n_surveys_aru,
        aru_independent_data,
        n_surveys_scores,
    ) in param_combinations:
        # Skip invalid combinations
        if not aru_independent_data and n_surveys_aru != n_surveys_scores:
            continue

        sim_params = SimulationParams(
            nsites=n_sites,  # type: ignore
            nsurveys_pc=n_surveys_pc,  # type: ignore
            nsurveys_aru=n_surveys_aru,  # type: ignore
            nsurveys_scores=n_surveys_scores,  # type: ignore
            beta0=beta0,  # type: ignore
            beta1=beta1,  # type: ignore
            p11=p11,  # type: ignore
            p_aru11=p_aru11,  # type: ignore
            p_aru01=p_aru01,  # type: ignore
            mu=mu,  # type: ignore
            sigma=sigma,  # type: ignore
            aru_data_independent_model=aru_independent_data,  # type: ignore
        )

        sim_args.append((sim_params, num_datasets, params.sim_name_for_data, name))

    print(
        f"Generating datasets for {len(sim_args)} parameter combinations across {num_processes} processes..."
    )

    # Parallel execution: Generate datasets in parallel
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = [
            executor.submit(_generate_datasets_for_params, args) for args in sim_args
        ]

        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Generating & saving datasets",
        ):
            sim_params, datasets = future.result()

            # Write results immediately as each parameter combination completes
            sim_param_id = db.insert_sim_params(study_id, sim_params)
            for data in datasets:
                db.insert_dataset(sim_param_id, data)
            db.commit()

    print("All datasets saved to database.")
    return study_id


def _run_single_simulation(
    args: Tuple[int, List[ModelNames], str],
) -> Tuple[int, List[Tuple[str, Any]]]:
    """
    Worker function: Run all models for a single dataset.

    Each process gets its own DB connection for reading only.
    Results are returned to main process for thread-safe writing.
    """
    dataset_id, model_names, _ = args

    # Each worker creates its own DB connection (read-only usage)
    db = SimulationsDB.create(DATABASE_PATH)
    dataset = db.get_dataset(dataset_id)
    params = db.get_sim_params(dataset.sim_param_id)

    results = []
    for model_name in model_names:
        # skip duplicates
        if db.get_run_ids_by_filters(
            params.study_id,
            runs_filter={"model_name": model_name, "dataset_id": dataset.id},
        ):
            continue
        # skip incompatible data/model combinations
        if (
            params.simulation_params.nsurveys_scores
            != params.simulation_params.nsurveys_aru
            and model_name == "single_year_jags_model_dependent"
        ):
            continue
        # we cannot run the independently generated aru/scores with the dependent model
        # because of the truncated normal distribution
        if (
            params.simulation_params.aru_data_independent_model
            and model_name == "single_year_jags_model_dependent"
        ):
            continue
        model = model_zoo.get_model_by_name(model_name)
        result = model.run_model(dataset.data)
        results.append((model_name, result))

    return dataset_id, results


def simulate_data(db: SimulationsDB, study_id: int, num_processes: int | None = None):
    """
    Run simulations in parallel and write results to database.

    Thread-safe: Each worker process runs simulations independently,
    then all results are written serially in the main process.

    Args:
        db: Database connection
        study_id: Study ID to run simulations for
        num_processes: Number of parallel processes (defaults to CPU count)
    """
    if num_processes is None:
        num_processes = cpu_count()

    dataset_ids = db.get_all_dataset_ids(study_id)
    study = db.get_study(study_id)

    print(f"Running {len(dataset_ids)} simulations across {num_processes} processes...")

    sim_args = [
        (dataset_id, study.study_params.models, study.study_params.sim_name_for_data)
        for dataset_id in dataset_ids
    ]

    # Parallel execution: ProcessPoolExecutor allows PyMC to use multiprocessing
    # Write results to DB as each simulation completes to minimize memory usage
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = [executor.submit(_run_single_simulation, args) for args in sim_args]

        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Running & saving simulations",
        ):
            dataset_id, model_results = future.result()

            # Write results immediately as each simulation completes
            for model_name, result in model_results:
                db.insert_run(dataset_id, model_name, result)
            db.commit()

    print("All results saved to database.")


if __name__ == "__main__":
    db = SimulationsDB.create(DATABASE_PATH)

    params = StudyParams(
        models=[
            "single_year_jags_model_dependent",
            "single_year_jags_model_independent",
        ],
        beta0=[-2.5, -1, 1, 2.5],
        beta1=[-1, 0.5],
        p11=[0.1, 0.5, 0.9],
        p_aru11=[0.1, 0.5, 0.9],
        p_aru01=[0, 0.05],
        mu=[(-2, -1.75), (-2, 0)],
        sigma=[(0.25, 1)],
        n_sites=[80, 200],
        n_surveys_pc=[3],
        n_surveys_aru=[24],
        n_surveys_scores=[8, 24],
        sim_name_for_data="single_year_jags_model_dependent",
        aru_scores_independent_data=[True, False],
        threshold=[-1, 0, 1],
    )

    # study_id = make_datasets(
    #     db, params, "full_comb_params", num_datasets=100, num_processes=32
    # )
    # print(f"Created {len(db.get_all_sim_param_ids(study_id))} parameter combinations")

    # Run simulations in parallel and save to database (thread-safe)
    study_id = 1
    simulate_data(db, study_id, num_processes=48)
