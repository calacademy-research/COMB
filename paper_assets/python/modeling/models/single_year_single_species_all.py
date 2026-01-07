from arviz import InferenceData
from caples_data import COMBData

from .model_iterface import CombinedModelInterface, SimulationParams, standardize
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
        burn_norm = standardize(data.covariates["caples"])
        # this is a single year, single species model so we need to extract the correct dimensions
        burn = burn_norm
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
        with pm.Model():
            # Set seed if provided
            if sim_params.seed is not None:
                np.random.seed(sim_params.seed)

            # ---------------------
            # Covariate
            # ---------------------
            covar = pm.Normal("covar", mu=0, sigma=1, shape=sim_params.nsites)

            # ---------------------
            # Occupancy
            # ---------------------
            logit_psi = sim_params.beta0 + sim_params.beta1 * covar
            psi = pm.Deterministic("psi", pm.math.sigmoid(logit_psi))
            z = pm.Bernoulli("z", p=psi, shape=sim_params.nsites)

            # ---------------------
            # Point count (PC)
            # ---------------------
            p_pc = z * sim_params.p11
            pm.Bernoulli(
                "y_ind",
                p=p_pc[:, None],
                shape=(sim_params.nsites, sim_params.nsurveys_pc),
            )

            # ---------------------
            # ARU detections (independent)
            # ---------------------
            p_aru = z * sim_params.p_aru11 + (1 - z) * sim_params.p_aru01
            pm.Bernoulli(
                "y_aru",
                p=p_aru[:, None],
                shape=(sim_params.nsites, sim_params.nsurveys_aru),
            )

            # ---------------------
            # Gaussian mixture scores
            # ---------------------
            mu_score = pm.math.switch(z, sim_params.mu[1], sim_params.mu[0])
            sigma_score = pm.math.switch(z, sim_params.sigma[1], sim_params.sigma[0])

            pm.Normal(
                "scores",
                mu=mu_score[:, None],  # type: ignore
                sigma=sigma_score[:, None],  # type: ignore
                shape=(sim_params.nsites, sim_params.nsurveys_scores),
            )

            # Sample from prior predictive
            prior_pred = pm.sample_prior_predictive(
                samples=1, random_seed=sim_params.seed
            )

        # Extract the simulated data (take the first sample from chain 0)
        covar_sim = prior_pred["prior"]["covar"].values[0, 0]
        y_ind_sim = prior_pred["prior"]["y_ind"].values[0, 0]
        y_aru_sim = prior_pred["prior"]["y_aru"].values[0, 0]
        scores_sim = prior_pred["prior"]["scores"].values[0, 0]

        return COMBData(
            n_sites=sim_params.nsites,
            n_years=1,
            n_species=1,
            y_aru=y_aru_sim.reshape(1, 1, y_aru_sim.shape[0], y_aru_sim.shape[1]),
            y_index=y_ind_sim.reshape(1, 1, y_ind_sim.shape[0], y_ind_sim.shape[1]),
            y_pc=np.zeros(0),
            n_surveys_pc=sim_params.nsurveys_pc,
            n_surveys_aru=sim_params.nsurveys_aru,
            covariates={"caples": covar_sim},
            date_pc=np.zeros(0),
            date_aru=np.zeros(0),
            time_aru=np.zeros(0),
            time_pc=np.zeros(0),
            scores=scores_sim.reshape(1, 1, scores_sim.shape[0], scores_sim.shape[1]),
        )
