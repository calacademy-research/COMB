import abc
from caples_data import COMBData
from arviz import InferenceData
from dataclasses import dataclass
import numpy as np


def normalize(x: np.ndarray) -> np.ndarray:
    return (x - np.nanmean(x)) / np.nanstd(x)


@dataclass
class SimulationParams:
    nsites: int = 50
    nsurveys_pc: int = 3
    nsurveys_aru: int = 5
    nsurveys_scores: int = 5
    beta0: float = -0.5
    beta1: float = 1.0
    p11: float = 0.6
    p_aru11: float = 0.7
    p_aru01: float = 0.1
    mu: tuple = (-1.0, 2.0)
    sigma: tuple = (1.0, 1.0)
    seed: int | None = None


class CombinedModelInterface(abc.ABC):
    @classmethod
    @abc.abstractmethod
    def run_model(cls, data: COMBData) -> InferenceData:
        """
        The main method to run the given model with COMBData. Returns an InferenceData object containing the samples.
        """

    @classmethod
    @abc.abstractmethod
    def simulate_data(cls, sim_params: SimulationParams) -> COMBData:
        """
        Simulate date for the model with the given params.
        """
