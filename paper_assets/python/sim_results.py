import arviz as az
import pandas as pd
import numpy as np
from sim_data import SimParams
from typing import List, Union
from pathlib import Path
import json


class SimResults:
    """
    SimResults class organizes the results from multiple simulations
    """

    def __init__(self, sim_params: SimParams):
        """
        Initializes a SimResults object given the simulation parameters.

        Args:
            sim_params (SimParams): The simulation parameters.
        """

        self.sim_params = sim_params
        # samples_dict is a list of dictionaries,
        # where each dictionary contains the samples from a single simulation
        self.samples_list: List[dict] = []

    def append_samples(self, samples: dict):
        """
        Append samples from a single simulation to the samples_dict

        :param samples: dictionary of samples from a single simulation
        """
        self.samples_list.append(samples)

    def get_summary(self) -> dict:
        """
        Get the summary for all of the samples

        :return summary: A dictionary of DataFrames, where keys are parameter names and
        each DataFrame contains the summary statistics for a single parameter
        """
        summary_dict = {}

        for samples in self.samples_list:
            self._summarize_single_simulation(samples, summary_dict)

        return summary_dict

    def _summarize_single_simulation(
        self,
        samples: dict,
        samples_dict: dict,
    ) -> None:
        """
        Summarize the samples from a single simulation

        :param samples: dictionary of samples from a single simulation
        :param samples_dict: dictionary of DataFrames, where keys are parameter names and
        each DataFrame contains the summary statistics for a single parameter
        """
        pyjags_data = az.from_pyjags(posterior=samples)
        summary = az.summary(pyjags_data)
        for param_name, single_param_summary in summary.iterrows():
            if param_name not in samples_dict:
                samples_dict[param_name] = pd.DataFrame()
            single_param_summary_df = single_param_summary.to_frame().T
            samples_dict[param_name] = pd.concat(
                [samples_dict[param_name], single_param_summary_df], ignore_index=True
            )

    def save_summary(self, directory: Union[Path, str]):
        """
        Save the summary statistics to a directory

        :param directory: The directory to save the summary statistics
        """
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        metadata = self.sim_params.__dict__
        for key, value in metadata.items():
            if isinstance(value, np.ndarray):
                metadata[key] = value.tolist()

        with open(directory / "metadata.json", "w") as f:
            f.write(json.dumps(metadata, indent=4))

        summary_dict = self.get_summary()

        for param_name, summary_df in summary_dict.items():
            summary_df.to_csv(directory / f"{param_name}_summary.csv")
