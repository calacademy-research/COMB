from arviz import InferenceData
from caples_data import COMBData

from .model_iterface import CombinedModelInterface, SimulationParams, normalize
import numpy as np
import pymc as pm


class SingleYearSingleSpeciesNoScoresNoPC(CombinedModelInterface):
    @classmethod
    def run_model(cls, data: COMBData) -> InferenceData:
        burn_norm = normalize(data.covariates["caples"])
        # this is a single year, single species model so we need to extract the correct dimensions
        burn = burn_norm[0]
        y_aru = data.y_aru[0, 0]

        print("burn", burn.shape)
        print("y_aru", y_aru.shape)

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
                    p_aru11,
                    p_aru01,
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
            # ARU detections (independent)
            # ---------------------
            p_aru = z * sim_params.p_aru11 + (1 - z) * sim_params.p_aru01
            pm.Bernoulli(
                "y_aru",
                p=p_aru[:, None],
                shape=(sim_params.nsites, sim_params.nsurveys_aru),
            )

            # Sample from prior predictive
            prior_pred = pm.sample_prior_predictive(
                samples=1, random_seed=sim_params.seed
            )

        # Extract the simulated data (take the first sample from chain 0)
        covar_sim = prior_pred["prior"]["covar"].values[0, 0]
        y_aru_sim = prior_pred["prior"]["y_aru"].values[0, 0]

        return COMBData(
            n_sites=sim_params.nsites,
            n_years=1,
            n_species=1,
            y_aru=y_aru_sim.reshape(1, 1, y_aru_sim.shape[0], y_aru_sim.shape[1]),
            y_index=np.zeros(0),
            y_pc=np.zeros(0),
            n_surveys_pc=sim_params.nsurveys_pc,
            n_surveys_aru=sim_params.nsurveys_aru,
            covariates={"burn": covar_sim.reshape(1, -1)},
            date_pc=np.zeros(0),
            date_aru=np.zeros(0),
            time_aru=np.zeros(0),
            time_pc=np.zeros(0),
            scores=np.zeros(0),
        )
