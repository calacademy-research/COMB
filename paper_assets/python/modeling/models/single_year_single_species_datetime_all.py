from arviz import InferenceData
from caples_data import COMBData
import numpy as np

from .model_iterface import CombinedModelInterface, normalize
import pymc as pm


class SingleYearSingleSpeciesAllDateTime(CombinedModelInterface):
    @classmethod
    def run_model(cls, data: COMBData) -> InferenceData:
        burn_norm = normalize(data.covariates["caples"])
        # this is a single year, single species model so we need to extract the correct dimensions
        burn = burn_norm[0]
        y_ind = data.y_index[0, 0]
        y_aru = data.y_aru[0, 0]
        scores = data.scores[0, 0]

        date_pc = normalize(data.date_pc[0, 0])
        date_aru = normalize(data.date_aru[0, 0])
        time_pc = normalize(data.time_pc[0, 0])
        time_aru = normalize(data.time_aru[0, 0])

        # Replace NaN values with 0 (normalized mean) to prevent NaN propagation
        date_pc = np.nan_to_num(date_pc, nan=0.0)
        time_pc = np.nan_to_num(time_pc, nan=0.0)
        date_aru = np.nan_to_num(date_aru, nan=0.0)
        time_aru = np.nan_to_num(time_aru, nan=0.0)

        print("burn", burn.shape)
        print("y_ind", y_ind.shape)
        print("y_aru", y_aru.shape)
        print("scores", scores.shape)
        print("date_pc", date_pc.shape)
        print("date_aru", date_aru.shape)
        print("time_pc", time_pc.shape)
        print("time_aru", time_aru.shape)

        with pm.Model() as model:
            # ---------------------
            # Occupancy
            # ---------------------
            beta0 = pm.Normal("beta0", mu=0, sigma=2)
            beta1 = pm.Normal("beta1", mu=0, sigma=2)

            logit_psi = beta0 + beta1 * burn
            psi = pm.Deterministic("psi", pm.math.sigmoid(logit_psi))

            z = pm.Bernoulli("z", p=psi, shape=data.n_sites)

            # ---------------------
            # Point count (PC)
            # ---------------------
            alpha0 = pm.Uniform("alpha0", lower=-5, upper=5)
            alpha1 = pm.Uniform("alpha1", lower=-5, upper=5)
            alpha2 = pm.Uniform("alpha2", lower=-5, upper=5)
            alpha3 = pm.Uniform("alpha3", lower=-5, upper=5)

            p11 = pm.Deterministic(
                "p11",
                pm.math.sigmoid(
                    alpha0
                    + alpha1 * time_pc
                    + alpha2 * time_pc * time_pc
                    + alpha3 * date_pc
                ),
            )

            p_pc = z[:, None] * p11
            pm.Bernoulli(
                "y_pc",
                p=p_pc,
                observed=y_ind,
            )

            # ---------------------
            # ARU detections
            # ---------------------
            p_aru11 = pm.Beta("p_aru11", alpha=1, beta=1)
            p_aru01 = pm.Beta("p_aru01", alpha=1, beta=1)

            p_aru = z[:, None] * p_aru11 + (1 - z[:, None]) * p_aru01
            pm.Bernoulli(
                "y_aru",
                p=p_aru,
                observed=y_aru,
            )

            # ---------------------
            # Gaussian mixture scores
            # ---------------------
            mu = pm.Normal("mu", mu=[-1, 1], sigma=2, shape=2)
            sigma = pm.HalfNormal("sigma", sigma=1, shape=2)

            mu_score = pm.math.switch(z, mu[1], mu[0])
            sigma_score = pm.math.switch(z, sigma[1], sigma[0])

            pm.Normal(
                "scores",
                mu=mu_score[:, None],  # type: ignore
                sigma=sigma_score[:, None],  # type: ignore
                observed=scores,
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
            model.debug()

            step_z = pm.BinaryGibbsMetropolis(vars=[z])
            step_cont = pm.NUTS(
                vars=[
                    beta0,
                    beta1,
                    alpha0,
                    alpha1,
                    alpha2,
                    alpha3,
                    p_aru11,
                    p_aru01,
                    mu,
                    sigma,
                ],
                target_accept=0.9,
            )

            trace = pm.sample(
                step=[step_z, step_cont],
                progressbar=True,
            )

        return trace
