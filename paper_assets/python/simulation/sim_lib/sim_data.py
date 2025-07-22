from dataclasses import dataclass, field
from typing import Tuple
import numpy as np
import hashlib
import json


@dataclass
class ModelParams:
    """
    ModelParams class hold the parameters that control the model, such as
    whether to include the covar model, scores model, etc...
    """

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
    aru_scores_independent_model: bool = True  # do not change


@dataclass
class DataParams:
    """
    DataParams class holds the parameters that control the generated data,
    such as psi, n_surveys, etc...
    """

    psi: float = 0.5  # only used if include_covar_model is False
    p11: float = 0.5
    p_aru11: float = 0.5
    p_aru01: float = 0.05
    beta0: float = 0  # only used if include_covar_model is True
    beta1: float = 0  # only used if include_covar_model is True
    mu: Tuple[float, float] = (-2, 2)
    sigma: Tuple[float, float] = (0.5, 2)
    nsites: int = 100
    # these two below must be the same if aru_scores_independent_model is False
    nsurveys_aru: int = 24
    nsurveys_scores: int = 24
    nsurveys_pc: int = 3
    covar_continuous: bool = True
    covar_prob: float = 0.5
    threshold: float = 0.5

    tau: Tuple[float, ...] = (4, 0.25)
    siteid: Tuple[int, ...] = (1, 2, 3)
    nsamples: int = -1
    covar: np.ndarray = field(default_factory=lambda: np.array([]))

    include_aru_data: bool = True
    include_pc_data: bool = True
    include_scores_data: bool = True
    include_covar_data: bool = True

    def __post_init__(self):
        # We can only have a diff number of scores and aru samples if the aru model is independent
        # if not self.aru_scores_independent_model:
        #     assert self.nsurveys_aru == self.nsurveys_scores, (
        #         "nsurveys_aru and nsurveys_scores must be the same if aru_scores_independent_model is False"
        #     )
        self.tau = tuple([1 / (s * s) for s in self.sigma])
        self.siteid = tuple(list(range(1, self.nsites + 1)))
        self.nsamples = len(self.siteid)


@dataclass
class SimParams:
    """
    SimParams class to hold simulation parameters

    These are the parameters that will be used to simulate data and run the model.

    """

    model_params: ModelParams
    data_params: DataParams

    def keys(self):
        return self.__dict__.keys()

    def __hash__(self):
        """
        Hashes the parameters using a frozen set of the parameters
        """
        data_dict = self.data_params.__dict__
        model_dict = self.model_params.__dict__
        d = {"data": data_dict, "model": model_dict}
        h = hashlib.md5((json.dumps(d, sort_keys=True, default=str)).encode()).digest()

        return int.from_bytes(h[:8], "big")
