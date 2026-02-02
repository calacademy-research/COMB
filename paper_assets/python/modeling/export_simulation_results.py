import polars as pl
from modeling.simulation.simulations_db import SimulationsDB
import arviz as az
import numpy as np
from tqdm import tqdm
from collections import defaultdict


"""
Script to export the simulation results into a summary format used in the paper
"""

schema = pl.Schema(
    {
        "run_id": pl.Int32,
        "dataset_id": pl.Int32,
        "params_id": pl.Int32,
        "mean": pl.Float32,
        "median": pl.Float32,
        "sd": pl.Float32,
        "hdi_3%": pl.Float32,
        "hdi_97%": pl.Float32,
        "mcse_mean": pl.Float32,
        "mcse_sd": pl.Float32,
        "ess_bulk": pl.Float32,
        "ess_tail": pl.Float32,
        "r_hat": pl.Float32,
        "beta0": pl.Float32,
        "beta1": pl.Float32,
        "p11": pl.Float32,
        "p_aru11": pl.Float32,
        "p_aru01": pl.Float32,
        "mu0": pl.Float32,
        "mu1": pl.Float32,
        "sigma0": pl.Float32,
        "sigma1": pl.Float32,
        "nsites": pl.Int32,
        "nsurveys_pc": pl.Int32,
        "nsurveys_aru": pl.Int32,
        "nsurveys_scores": pl.Int32,
        "threshold": pl.Float32,
        "aru_model_independent": pl.Boolean,
        "aru_data_independent": pl.Boolean,
    }
)


def median(x):
    return np.median(x)


if __name__ == "__main__":
    study_id = 1
    db = SimulationsDB.create("data/comb_simulations")

    run_ids = db.get_run_ids_by_filters(study_id)

    # each param has its own df, so the param (beta0, p_aru11, ...)
    # gets mapped to a list of rows
    summary_rows: dict[str, list] = defaultdict(list)

    for run_id in tqdm(run_ids[:5]):
        run = db.get_run(run_id)
        dataset_id = run.dataset_id
        dataset = db.get_dataset(dataset_id)
        params = db.get_sim_params(dataset.sim_param_id)

        result_summary: pl.DataFrame = pl.from_pandas(
            az.summary(run.results, stat_funcs={"median": median}),
            include_index=True,
        ).rename({"None": "param"})

        is_independent_model = run.model_name == "single_year_jags_model_independent"

        for row in result_summary.iter_rows(named=True):
            param_name = row["param"]

            summary_rows[param_name].append(
                {
                    "run_id": run_id,
                    "dataset_id": dataset_id,
                    "params_id": dataset.sim_param_id,
                    "mean": row["mean"],
                    "median": row["median"],
                    "sd": row["sd"],
                    "hdi_3%": row["hdi_3%"],
                    "hdi_97%": row["hdi_97%"],
                    "mcse_mean": row["mcse_mean"],
                    "mcse_sd": row["mcse_sd"],
                    "ess_bulk": row["ess_bulk"],
                    "ess_tail": row["ess_tail"],
                    "r_hat": row["r_hat"],
                    "beta0": params.simulation_params.beta0,
                    "beta1": params.simulation_params.beta1,
                    "p11": params.simulation_params.p11,
                    "p_aru11": params.simulation_params.p_aru11,
                    "p_aru01": params.simulation_params.p_aru01,
                    "mu0": params.simulation_params.mu[0],
                    "mu1": params.simulation_params.mu[1],
                    "sigma0": params.simulation_params.sigma[0],
                    "sigma1": params.simulation_params.sigma[1],
                    "nsites": params.simulation_params.nsites,
                    "nsurveys_pc": params.simulation_params.nsurveys_pc,
                    "nsurveys_aru": params.simulation_params.nsurveys_aru,
                    "nsurveys_scores": params.simulation_params.nsurveys_scores,
                    "threshold": params.simulation_params.threshold,
                    "aru_model_independent": is_independent_model,
                    "aru_data_independent": params.simulation_params.aru_data_independent_model,
                }
            )

    for param, rows in summary_rows.items():
        df = pl.DataFrame(rows, schema=schema)

        df.write_csv(f"data/comb_simulations_outputs/{param}.csv")
