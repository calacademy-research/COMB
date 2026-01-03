import pymc as pm
import numpy as np
from caples_data.combine_aru_pc import COMBData, CombinedData, CombinedParams
import polars as pl

# Data (shapes matter)
# nsites
# nsurveys_pc
# nsurveys_aru
# nsurveys_scores

# covar:        (nsites,)
# y_pc:         (nsites, nsurveys_pc)
# y_aru:        (nsites, nsurveys_aru)
# scores:       (nsites, nsurveys_scores)
# siteid:       (nsites,)  # 0-based indexing


class CombinedModel:
    def __init__(self, data: COMBData):
        self.data = data
        # normalize covar data
        self.data.covariates["burn"] = (
            self.data.covariates["burn"] - np.mean(self.data.covariates["burn"])
        ) / np.std(self.data.covariates["burn"])

    def run_model(self):
        siteid = np.arange(self.data.n_sites)
        burn = self.data.covariates["burn"]

        with pm.Model():
            # ---------------------
            # Occupancy
            # ---------------------
            beta0 = pm.Normal("beta0", mu=0, sigma=2)
            beta1 = pm.Normal("beta1", mu=0, sigma=2)

            logit_psi = beta0 + beta1 * burn
            psi = pm.Deterministic("psi", pm.math.sigmoid(logit_psi))

            z = pm.Bernoulli("z", p=psi, shape=self.data.n_sites)

            # ---------------------
            # Point count (PC)
            # ---------------------
            p11 = pm.Beta("p11", alpha=1, beta=1)

            p_pc = z * p11
            pm.Bernoulli(
                "y_pc",
                p=p_pc[:, None],
                observed=self.data.y_pc,
            )

            # ---------------------
            # ARU detections
            # ---------------------
            p_aru11 = pm.Beta("p_aru11", alpha=1, beta=1)
            p_aru01 = pm.Beta("p_aru01", alpha=1, beta=1)

            p_aru = z * p_aru11 + (1 - z) * p_aru01
            pm.Bernoulli(
                "y_aru",
                p=p_aru[:, None],
                observed=self.data.y_aru,
            )

            # ---------------------
            # Gaussian mixture scores
            # ---------------------
            mu = pm.Normal("mu", mu=0, sigma=5, shape=2)
            sigma = pm.HalfNormal("sigma", sigma=2, shape=2)

            mu_score = pm.math.switch(z[siteid], mu[1], mu[0])
            sigma_score = pm.math.switch(z[siteid], sigma[1], sigma[0])

            pm.Normal(
                "scores",
                mu=mu_score[:, None],  # type: ignore
                sigma=sigma_score[:, None],  # type: ignore
                observed=self.data.scores,
            )

            # ---------------------
            # Derived quantities
            # ---------------------
            pm.Deterministic("mean_psi", psi.mean())
            pm.Deterministic("NOcc", z.sum())
            pm.Deterministic("PropOcc", z.mean())

            # ---------------------
            # SAMPLING (IMPORTANT)
            # ---------------------
            step_z = pm.Metropolis(vars=[z])
            step_cont = pm.NUTS(
                vars=[beta0, beta1, p11, p_aru11, p_aru01, mu, sigma],
                target_accept=0.9,
            )

            trace = pm.sample(
                step=[step_z, step_cont],
                progressbar=True,
            )

        return trace


if __name__ == "__main__":
    aru = pl.read_parquet(
        "/Users/mschulist/github/COMB/paper_assets/python/caples_data/data/outputs_agg_20251224_180943.parquet"
    ).filter(pl.col("point") != 0)

    pc = pl.read_csv(
        "/Users/mschulist/github/COMB/paper_assets/python/caples_data/data/PC_delinted.csv"
    ).with_columns(
        visit=pl.col("visit") - 1,
        DateTime=pl.col("DateTime").str.to_datetime("%Y-%m-%dT%H:%M:%SZ"),
    )

    spatial = pl.read_csv(
        "/Users/mschulist/github/COMB/paper_assets/python/spatial/data/burn_data_by_point.csv"
    )

    combined_params = CombinedParams(
        aru_species_col="label",
        aru_visit_limit=24,
        years=[2021],
        pc_species_col="birdCode_fk",
        pc_count_col="abun",
        pc_datetime_col="DateTime",
        pc_point_col="point_ID_fk",
        pc_visit_index_col="visit",
        species=["gockin"],
    )
    combined = CombinedData(
        aru,
        pc,
        spatial,
        combined_params,
    )

    model = CombinedModel(combined.combined_data)
    trace = model.run_model()

    trace.to_netcdf("data/trace.nc")
