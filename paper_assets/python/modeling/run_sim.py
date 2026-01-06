from simulation.study_params import StudyParams
from models import model_zoo
from models.model_iterface import SimulationParams
from models.model_zoo import ModelNames
from simulation.simulations_db import SimulationsDB
from datetime import datetime as dt
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from typing import Any, Tuple, List
from tqdm import tqdm


def make_datasets(
    db: SimulationsDB, params: StudyParams, name: str, num_datasets: int
) -> int:
    """Create study and generate datasets for all parameter combinations."""
    study_id = db.insert_study(name, dt.now(), params)
    model = model_zoo.get_model_by_name(params.sim_name_for_data)

    for beta0 in params.beta0:
        for beta1 in params.beta1:
            for p11 in params.p11:
                for p_aru11 in params.p_aru11:
                    for p_aru01 in params.p_aru01:
                        for mu in params.mu:
                            for sigma in params.sigma:
                                for n_sites in params.n_sites:
                                    for n_surveys_pc in params.n_surveys_pc:
                                        for n_surveys_aru in params.n_surveys_aru:
                                            for (
                                                n_surveys_scores
                                            ) in params.n_surveys_scores:
                                                sim_params = SimulationParams(
                                                    nsites=n_sites,
                                                    nsurveys_pc=n_surveys_pc,
                                                    nsurveys_aru=n_surveys_aru,
                                                    nsurveys_scores=n_surveys_scores,
                                                    beta0=beta0,
                                                    beta1=beta1,
                                                    p11=p11,
                                                    p_aru11=p_aru11,
                                                    p_aru01=p_aru01,
                                                    mu=mu,
                                                    sigma=sigma,
                                                )

                                                sim_param_id = db.insert_sim_params(
                                                    study_id, sim_params
                                                )

                                                for _ in range(num_datasets):
                                                    data = model.simulate_data(
                                                        sim_params
                                                    )
                                                    db.insert_dataset(
                                                        sim_param_id, data
                                                    )
    db.commit()
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
    db = SimulationsDB.create("data/simulations")
    dataset = db.get_dataset(dataset_id)

    results = []
    for model_name in model_names:
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
    db = SimulationsDB.create("data/simulations")

    params = StudyParams(
        models=[
            "single_year_single_species_all",
            "single_year_single_species_no_pc",
            "single_year_single_species_no_aru",
            "single_year_single_species_no_scores",
            "single_year_single_species_no_aru_no_pc",
            "single_year_single_species_no_scores_no_pc",
            "single_year_single_species_no_scores_no_aru",
        ],
        beta0=[-0.5],
        beta1=[-1, 0.5],
        p11=[0.1, 0.5],
        p_aru11=[0.1, 0.5],
        p_aru01=[0.05],
        mu=[(-2, -1.5), (-2, 0)],
        sigma=[(0.5, 1)],
        n_sites=[80, 200],
        n_surveys_pc=[3],
        n_surveys_aru=[24],
        n_surveys_scores=[24],
        sim_name_for_data="single_year_single_species_all",
    )

    study_id = make_datasets(db, params, "minimal_param_combs", num_datasets=1)
    print(f"Created {len(db.get_all_sim_param_ids(study_id))} parameter combinations")

    # Run simulations in parallel and save to database (thread-safe)
    simulate_data(db, study_id, num_processes=4)
