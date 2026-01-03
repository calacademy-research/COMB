from arviz import InferenceData
from caples_data import COMBData

from modeling.models.model_iterface import SimulationParams
from .model_iterface import CombinedModelInterface
import numpy as np
import pymc as pm

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


class SingleYearSingleSpeciesAll(CombinedModelInterface):
    @classmethod
    def run_model(cls, data: COMBData) -> InferenceData:
        burn = data.covariates["burn"] - np.mean(data.covariates["burn"]) / np.std(
            data.covariates["burn"]
        )
        # this is a single year, single species model so we need to extract the correct dimensions
        burn = data.covariates["burn"][0]
        y_ind = data.y_index[0, 0]
        y_aru = data.y_aru[0, 0]
        scores = data.scores[0, 0]

        print("burn", burn.shape)
        print("y_ind", y_ind.shape)
        print("y_aru", y_aru.shape)
        print("scores", scores.shape)

        with pm.Model():
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
            p11 = pm.Beta("p11", alpha=1, beta=1)

            p_pc = z * p11
            pm.Bernoulli(
                "y_pc",
                p=p_pc[:, None],
                observed=y_ind,
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
            step_z = pm.BinaryGibbsMetropolis(vars=[z])
            step_cont = pm.NUTS(
                vars=[
                    beta0,
                    beta1,
                    p11,
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

    @classmethod
    def simulate_data(cls, sim_params: SimulationParams) -> COMBData:
        if sim_params.seed is not None:
            np.random.seed(sim_params.seed)

        # ---------------------
        # Covariate
        # ---------------------
        covar = np.random.normal(0, 1, sim_params.nsites)

        # ---------------------
        # Occupancy
        # ---------------------
        logit_psi = sim_params.beta0 + sim_params.beta1 * covar
        psi = 1 / (1 + np.exp(-logit_psi))
        z = np.random.binomial(1, psi, size=sim_params.nsites)

        # ---------------------
        # Point count (PC)
        # ---------------------
        p_pc = z * sim_params.p11
        y_ind = np.random.binomial(
            1,
            p_pc[:, None],
            size=(sim_params.nsites, sim_params.nsurveys_pc),
        )

        # ---------------------
        # ARU detections (independent)
        # ---------------------
        p_aru = z * sim_params.p_aru11 + (1 - z) * sim_params.p_aru01
        y_aru = np.random.binomial(
            1,
            p_aru[:, None],
            size=(sim_params.nsites, sim_params.nsurveys_aru),
        )

        # ---------------------
        # Gaussian mixture scores
        # ---------------------
        scores = np.zeros((sim_params.nsites, sim_params.nsurveys_scores))
        for i in range(sim_params.nsites):
            k = z[i]  # 0 or 1
            scores[i, :] = np.random.normal(
                loc=sim_params.mu[k],
                scale=sim_params.sigma[k],
                size=sim_params.nsurveys_scores,
            )

        return COMBData(
            n_sites=sim_params.nsites,
            n_years=1,
            n_species=1,
            y_aru=y_aru,
            y_index=y_ind,
            y_pc=np.zeros(0),
            n_surveys_pc=sim_params.nsurveys_pc,
            n_surveys_aru=sim_params.nsurveys_aru,
            covariates={"burn": covar},
            date_pc=np.zeros(0),
            date_aru=np.zeros(0),
            time_aru=np.zeros(0),
            time_pc=np.zeros(0),
            scores=scores,
        )
