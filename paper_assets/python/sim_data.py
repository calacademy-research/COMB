from dataclasses import dataclass
from typing import Tuple
import numpy as np
from scipy.special import expit, logit
import scipy.stats as stats


@dataclass
class SimParams:
    p11: float = 0.5
    p_aru11: float = 0.5
    p_aru01: float = 0.05
    beta0: float = 0
    beta1: float = 0
    mu: Tuple[float, float] = (-2, 2)
    sigma: Tuple[float, float] = (0.5, 2)
    nsites: int = 100
    nsurveys_aru: int = 24
    nsurveys_pc: int = 3
    covar_continuous: bool = False
    covar_prob: float = 0.5
    threshold: float = 0.5

    tau: Tuple[float, ...] = (4, 0.25)
    siteid: Tuple[int, ...] = (1, 2, 3)
    nsamples: int = -1
    covar: np.ndarray = np.ndarray([])

    psi_prior: str = "dbeta(2,2)"
    beta0_prior: str = "dnorm(0, 0.5)"
    beta1_prior: str = "dnorm(0, 0.5)"
    p11_prior: str = "dbeta(2, 2)"
    mu1_prior: str = "dnorm(-2, 0.2)"
    mu2_prior: str = "dnorm(-2, 0.2)"
    sigma1_prior: str = "dunif(0.1, 5)"
    sigma2_prior: str = "dunif(0.1, 5)"
    p_aru11_prior: str = "dbeta(2, 2)"
    p_aru01_prior: str = "dbeta(1, 3)I(0, 1 - p_aru11)"
    na: int = 1000  # number of iterations to adapt
    ni: int = 8000  # number of iterations
    nt: int = 1  # thinning rate
    nb: int = 1000  # number of burn-in iterations
    nc: int = 6  # number of chains
    parallel: bool = True  # whether to run each chain in parallel

    include_aru_model: bool = True
    include_pc_model: bool = True
    include_scores_model: bool = True
    include_covar_model: bool = True
    aru_scores_independent_model: bool = True
    
    def keys(self):
        return self.__dict__.keys()

    def __post_init__(self):
        self.nsurveys_scores = self.nsurveys_aru
        self.tau = tuple([1 / (s * s) for s in self.sigma])
        self.siteid = tuple(list(range(1, self.nsites + 1)))
        self.nsamples = len(self.siteid)
        if self.covar_continuous:
            self.covar = np.random.normal(0, 1, self.nsites)
        else:
            self.covar = np.random.binomial(1, self.covar_prob, self.nsites)


# class SimData:
#     """
#     SimData class to generate simulated data
#     """

#     y_pc: np.ndarray = np.ndarray([])
#     y_aru: np.ndarray = np.ndarray([])
#     scores: np.ndarray = np.ndarray([])

#     def __init__(self, params: SimParams):
#         # params holds all of the simulation parameters
#         self.params = params

#         # below are the parameters that will be used when running the model
#         self.nsites = params.nsites
#         self.nsurveys_aru = params.nsurveys_aru
#         self.nsurveys_scores = params.nsurveys_aru
#         self.nsurveys_pc = params.nsurveys_pc
#         self.covar = params.covar
#         self.siteid = params.siteid

#     def gen_simulated_data(self):
#         # using covariates, we generate latent occupancy states
#         z = self._gen_latent_occupancy()
#         # generate PC data
#         self.y_pc = self._gen_pc_data(z)
#         self.y_aru = self._gen_aru_data(z)
#         self.scores = self._gen_aru_scores(self.y_aru)

#     def _gen_latent_occupancy(self):
#         """
#         Generate latent occupancy states
#         """
#         psi_c = expit(self.params.beta0 + self.params.beta1 * self.params.covar)
#         z = np.random.binomial(1, psi_c, self.params.nsites)
#         return z

#     def _gen_pc_data(self, z: np.ndarray) -> np.ndarray:
#         """
#         Generate point count data

#         Args:
#             z (numpy.ndarray): The latent occupancy states.

#         Returns:
#             y_pc (numpy.ndarray): The point count data.
#         """
#         y_pc = np.zeros((self.params.nsites, self.params.nsurveys_pc))

#         for i in range(self.params.nsites):
#             prob_y_eq_1 = self.params.p11 * z[i]
#             y_pc[i] = np.random.binomial(1, prob_y_eq_1, self.params.nsurveys_pc)

#         return y_pc

#     def _gen_aru_data(self, z: np.ndarray) -> np.ndarray:
#         """
#         Generate ARU data

#         Args:
#             z (numpy.ndarray): The latent occupancy states.

#         Returns:
#             y_aru (numpy.ndarray): The ARU data.
#         """
#         y_aru = np.zeros((self.params.nsites, self.params.nsurveys_aru))

#         for i in range(self.params.nsites):
#             p_aru = z[i] * self.params.p_aru11 + self.params.p_aru01
#             y_aru[i] = np.random.binomial(1, p_aru, self.params.nsurveys_aru)

#         return y_aru

#     def _gen_aru_scores(self, y_aru: np.ndarray) -> np.ndarray:
#         """
#         Generate ARU scores

#         Args:
#             y_aru (numpy.ndarray): The binary ARU data.

#         Returns:
#             y_aru_scores (numpy.ndarray): The ARU scores.
#         """

#         scores = np.zeros((self.params.nsites, self.params.nsurveys_aru))
#         for i in range(self.params.nsites):
#             for j in range(self.params.nsurveys_aru):
#                 if y_aru[i, j] == 0:
#                     min = -10
#                     max = self.params.threshold
#                     scores[i, j] = stats.truncnorm.rvs(
#                         min, max, loc=self.params.mu[0], scale=self.params.sigma[0]
#                     )
#                 elif y_aru[i, j] == 1:
#                     min = self.params.threshold
#                     max = 10
#                     scores[i, j] = stats.truncnorm.rvs(
#                         min, max, loc=self.params.mu[1], scale=self.params.sigma[1]
#                     )
#                 else:
#                     raise ValueError("Invalid y_aru value, must be 0 or 1")

#         return scores
