import numpy as np
from typing import Dict, Any


def simulate_combined_data(
    nsites: int = 50,
    nsurveys_pc: int = 3,
    nsurveys_aru: int = 5,
    nsurveys_scores: int = 5,
    beta0: float = -0.5,
    beta1: float = 1.0,
    p11: float = 0.6,
    p_aru11: float = 0.7,
    p_aru01: float = 0.1,
    mu: tuple = (-1.0, 2.0),
    sigma: tuple = (1.0, 1.0),
    seed: int | None = None,
) -> Dict[str, Any]:
    if seed is not None:
        np.random.seed(seed)

    # ---------------------
    # Covariate
    # ---------------------
    covar = np.random.normal(0, 1, nsites)

    # ---------------------
    # Occupancy
    # ---------------------
    logit_psi = beta0 + beta1 * covar
    psi = 1 / (1 + np.exp(-logit_psi))
    z = np.random.binomial(1, psi, size=nsites)

    # ---------------------
    # Point count (PC)
    # ---------------------
    p_pc = z * p11
    y_pc = np.random.binomial(
        1,
        p_pc[:, None],
        size=(nsites, nsurveys_pc),
    )

    # ---------------------
    # ARU detections (independent)
    # ---------------------
    p_aru = z * p_aru11 + (1 - z) * p_aru01
    y_aru = np.random.binomial(
        1,
        p_aru[:, None],
        size=(nsites, nsurveys_aru),
    )

    # ---------------------
    # Gaussian mixture scores
    # ---------------------
    scores = np.zeros((nsites, nsurveys_scores))
    for i in range(nsites):
        k = z[i]  # 0 or 1
        scores[i, :] = np.random.normal(
            loc=mu[k],
            scale=sigma[k],
            size=nsurveys_scores,
        )

    return {
        "nsites": nsites,
        "covar": covar,
        "y_pc": y_pc,
        "y_aru": y_aru,
        "scores": scores,
    }
