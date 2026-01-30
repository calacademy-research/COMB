from arviz import InferenceData
from modeling.caples_data import COMBData
import numpy as np

from .model_iterface import CombinedModelInterface, normalize
import pymc as pm


class MultiYearMultiSpeciesAllMarg(CombinedModelInterface):
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

        S, T, I, J = y_aru.shape

        with pm.Model() as model:
            # ======================================================
            # INITIAL OCCUPANCY
            # ======================================================
            beta0 = pm.Normal("beta0", 0, 1, shape=S)
            psi0 = pm.Deterministic(
                "psi_0", pm.math.sigmoid(beta0[:, None]).repeat(I, axis=1)
            )
            psi_list = [psi0]

            # ======================================================
            # DYNAMICS
            # ======================================================
            phi0 = pm.Normal("phi0", 0, 1, shape=S)

            phi_caples_0 = pm.Normal("phi_caples_0", 0, 1, shape=S)
            phi_caples_1 = pm.Normal("phi_caples_1", 0, 1, shape=S)

            phi_caldor_0 = pm.Normal("phi_caldor_0", 0, 1, shape=S)
            phi_caldor_1 = pm.Normal("phi_caldor_1", 0, 1, shape=S)

            gamma0 = pm.Normal("gamma0", 0, 1, shape=S)
            gamma = pm.math.sigmoid(gamma0[:, None])

            for t in range(1, T):
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

                psi_prev = psi_list[-1]
                psi_t = psi_prev * phi + (1 - psi_prev) * gamma
                psi_list.append(psi_t)

            psi = pm.Deterministic("psi", pm.math.stack(psi_list, axis=1))  # (S,T,I)

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
            p_obs_pc = psi[..., None] * p_pc

            pm.Bernoulli("y_pc", p=p_obs_pc, observed=y_ind)

            # ======================================================
            # ARU DETECTION (FALSE POSITIVES)
            # ======================================================
            p_aru11 = pm.Beta("p_aru11", 1, 1, shape=S)
            p_aru01 = pm.Beta("p_aru01", 1, 1, shape=S)

            p_obs_aru = (
                psi[..., None] * p_aru11[:, None, None, None]
                + (1 - psi[..., None]) * p_aru01[:, None, None, None]
            )

            pm.Bernoulli("y_aru", p=p_obs_aru, observed=y_aru)

            # ======================================================
            # SCORE MIXTURE MODEL (MARGINALIZED)
            # ======================================================
            mu = pm.Normal("mu", [-1, 1], 2, shape=(S, 2))
            sigma = pm.HalfNormal("sigma", 1, shape=(S, 2))

            logp0 = pm.logp(
                pm.Normal.dist(
                    mu=mu[:, 0][:, None, None, None],
                    sigma=sigma[:, 0][:, None, None, None],
                ),
                scores,
            )

            logp1 = pm.logp(
                pm.Normal.dist(
                    mu=mu[:, 1][:, None, None, None],
                    sigma=sigma[:, 1][:, None, None, None],
                ),
                scores,
            )

            pm.Potential(
                "score_mixture",
                pm.math.logsumexp(
                    pm.math.stack(
                        [
                            pm.math.log(1 - psi[..., None]) + logp0,
                            pm.math.log(psi[..., None]) + logp1,
                        ],
                        axis=0,
                    ),
                    axis=0,
                ),
            )

            # ======================================================
            # DERIVED
            # ======================================================
            pm.Deterministic("mean_caples_resilience", phi_caples_1.mean())
            pm.Deterministic("mean_caldor_resilience", phi_caldor_1.mean())

            trace = pm.sample(
                target_accept=0.9,
                chains=4,
                cores=4,
            )

        return trace
