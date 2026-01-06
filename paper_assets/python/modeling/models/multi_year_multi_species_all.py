from arviz import InferenceData
from caples_data import COMBData
import numpy as np

from .model_iterface import CombinedModelInterface, normalize
import pymc as pm


class MultiYearMultiSpeciesAll(CombinedModelInterface):
    @classmethod
    def run_model(cls, data: COMBData) -> InferenceData:
        y_ind = data.y_index
        y_aru = data.y_aru
        scores = data.scores

        date_pc = np.nan_to_num(data.date_pc)
        time_pc = np.nan_to_num(data.time_pc)
        inside_caples = data.covariates["inside_caples"]
        inside_caldor = data.covariates["inside_caldor"]
        severity_caples = normalize(data.covariates["caples"])
        severity_caldor = normalize(data.covariates["caldor"])
        tsb_caples = data.covariates["time_since_caples"]
        tsb_caldor = data.covariates["time_since_caldor"]

        # species, time, points, visits
        S, T, I, J = y_aru.shape

        with pm.Model() as model:
            # ======================================================
            # INITIAL OCCUPANCY (pre-fire baseline)
            # ======================================================
            beta0 = pm.Normal("beta0", 0, 1, shape=S)
            psi0 = pm.math.sigmoid(beta0[:, None])
            z0 = pm.Bernoulli("z_0", psi0, shape=(S, I))
            z_list = [z0]

            # ======================================================
            # PERSISTENCE / RESILIENCE PARAMETERS
            # ======================================================
            # Baseline persistence (no fire)
            phi0 = pm.Normal("phi0", 0, 1, shape=S)

            # Immediate burn impact (first post-fire year)
            phi_caples_0 = pm.Normal("phi_caples_0", 0, 1, shape=S)
            phi_caldor_0 = pm.Normal("phi_caldor_0", 0, 1, shape=S)

            # Recovery slopes (per year since burn)
            phi_caples_1 = pm.Normal("phi_caples_1", 0, 1, shape=S)
            phi_caldor_1 = pm.Normal("phi_caldor_1", 0, 1, shape=S)

            # Colonization
            gamma0 = pm.Normal("gamma0", 0, 1, shape=S)

            # ======================================================
            # DYNAMIC OCCUPANCY
            # ======================================================
            for t in range(1, T):
                # Burn has an effect ONLY if time-since-burn > 0
                post_caples = pm.math.gt(tsb_caples[:, t], 0)
                post_caldor = pm.math.gt(tsb_caldor[:, t], 0)

                caples_effect = (
                    inside_caples[None, :]
                    * post_caples[None, :]
                    * severity_caples[None, :]
                    * (phi_caples_0[:, None] + phi_caples_1[:, None] * tsb_caples[:, t])
                )

                caldor_effect = (
                    inside_caldor[None, :]
                    * post_caldor[None, :]
                    * severity_caldor[None, :]
                    * (phi_caldor_0[:, None] + phi_caldor_1[:, None] * tsb_caldor[:, t])
                )

                logit_phi = phi0[:, None] + caples_effect + caldor_effect
                phi = pm.math.sigmoid(logit_phi)

                gamma = pm.math.sigmoid(gamma0[:, None])

                mu_z = z_list[t - 1] * phi + (1 - z_list[t - 1]) * gamma
                z_t = pm.Bernoulli(f"z_{t}", mu_z, shape=(S, I))
                z_list.append(z_t)

            z = pm.Deterministic("z", pm.math.stack(z_list, axis=1))

            # ======================================================
            # POINT COUNT DETECTION
            # ======================================================
            alpha0 = pm.Normal("alpha0", 0, 1, shape=S)
            alpha1 = pm.Normal("alpha1", 0, 1, shape=S)
            alpha2 = pm.Normal("alpha2", 0, 1, shape=S)
            alpha3 = pm.Normal("alpha3", 0, 1, shape=S)

            logit_p_pc = (
                alpha0[:, None, None, None]
                + alpha1[:, None, None, None] * time_pc
                + alpha2[:, None, None, None] * time_pc**2
                + alpha3[:, None, None, None] * date_pc
            )

            p_pc = pm.math.sigmoid(logit_p_pc)
            pm.Bernoulli("y_pc", p=z[..., None] * p_pc, observed=y_ind)

            # ======================================================
            # ARU DETECTION
            # ======================================================
            p_aru11 = pm.Beta("p_aru11", 1, 1, shape=S)
            p_aru01 = pm.Beta("p_aru01", 1, 1, shape=S)

            p_aru = (
                z[..., None] * p_aru11[:, None, None, None]
                + (1 - z[..., None]) * p_aru01[:, None, None, None]
            )

            pm.Bernoulli("y_aru", p=p_aru, observed=y_aru)

            # ======================================================
            # SCORE MIXTURE MODEL
            # ======================================================
            mu = pm.Normal("mu", [-1, 1], 2, shape=(S, 2))
            sigma = pm.HalfNormal("sigma", 1, shape=(S, 2))

            mu_score = pm.math.switch(z[..., None], mu[:, 1], mu[:, 0])
            sigma_score = pm.math.switch(z[..., None], sigma[:, 1], sigma[:, 0])

            pm.Normal("scores", mu=mu_score, sigma=sigma_score, observed=scores)

            # ======================================================
            # DERIVED ECOLOGICAL QUANTITIES
            # ======================================================
            pm.Deterministic("mean_caples_resilience", phi_caples_1.mean())
            pm.Deterministic("mean_caldor_resilience", phi_caldor_1.mean())

            # ======================================================
            # SAMPLING
            # ======================================================
            step_z = pm.BinaryGibbsMetropolis(vars=z_list)
            step_cont = pm.NUTS(target_accept=0.9)

            model.debug()

            trace = pm.sample(
                step=[step_z, step_cont],
                chains=16,
                cores=16,
                init="adapt_diag",
                progressbar=True,
            )

        return trace
