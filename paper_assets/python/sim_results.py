import arviz as az
import pandas as pd
import numpy as np
from sim_data import SimParams
from typing import List, Union
from pathlib import Path
import matplotlib.pyplot as plt
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

        self.summary: dict = {}

    def append_samples(self, samples: dict):
        """
        Append samples from a single simulation to the samples_dict

        Args:
            samples (dict): dictionary of samples from a single simulation
        """
        self.samples_list.append(samples)

    def get_summary(self) -> dict:
        """
        Get the summary for all of the samples

        Returns:
            dict: A dictionary of DataFrames, where keys are parameter names and each
            DataFrame contains the summary statistics for a single parameter
        """
        summary_dict = {}

        for samples in self.samples_list:
            self._summarize_single_simulation(samples, summary_dict)

        self.summary = summary_dict

        return summary_dict

    def _summarize_single_simulation(
        self,
        samples: dict,
        samples_dict: dict,
    ) -> None:
        """
        Summarize the samples from a single simulation

        Args:
            samples (dict): dictionary of samples from a single simulation
            samples_dict (dict): dictionary of DataFrames, where keys are parameter names and
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

        TODO: add ability to save all chains (compressed), not just the summary

        Args:
            directory (str | Path): The directory to save the summary statistics
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

    @staticmethod
    def load(directory: Union[str, Path]) -> "SimResults":
        """
        Load SimResults from directory

        You can save SimResults using the `save_summary` method

        Args:
            directory (str | Path)
        """
        # TODO: make this method
        pass
        return SimResults(SimParams())

    def plot_summary_hists(self, params: Union[List[str], None] = None):
        """
        Makes a histogram for each parameter of the mean value of that parameter

        Passing None to params will make a histogram for every parameter

        Args:
            params (List[str]): list of parameters to make histograms for
        """
        if not self.summary:
            self.get_summary()

        summary = self.summary

        figs = []
        for param_name, summary_df in summary.items():
            sim_param_value = self._get_sim_param_value(param_name)
            summary_values = summary_df["mean"]
            fig = plt.figure()
            plt.vlines(x=sim_param_value, ymin=0, ymax=max(summary_values))
            plt.hist(summary_values)
            plt.title(param_name)
            figs.append(fig)
        
        return figs

    def _get_sim_param_value(self, sim_param_name: str) -> Union[float, int]:
        """
        Helper method to get sim params values

        We cannot just query the SimParams with the `param_name` because
        params that are lists (such as mu and sigma) have the form
        `mu[0]` in the summary, but are lists (or tuples) in the SimParams
        class.

        Args:
            sim_param_name (str): name of the parameter as found in the summary

        Returns:
            (float | int): the value of the SimParam
        """
        sim_params = self.sim_params.__dict__
        if "[" not in sim_params and "]" not in sim_params and sim_param_name[-1] == "]":
            return sim_params[sim_param_name]
        name = "".join(sim_param_name.split("[")[0:-1])
        index = sim_param_name.split("[")[-1]
        try:
            return sim_params[name][index]
        except ValueError:
            raise ValueError(f"name {name} and index {index} not found")
