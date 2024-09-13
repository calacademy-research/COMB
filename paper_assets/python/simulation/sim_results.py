import math
import os
import arviz as az
import pandas as pd
import numpy as np
from sim_data import SimParams
from typing import List, Union
from pathlib import Path
import matplotlib.pyplot as plt
import json
import inspect


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
        Append samples from a single simulation to the samples_list

        Args:
            samples (dict): dictionary of samples from a single simulation
        """
        self.samples_list.append(samples)

    def get_summary(self, reindex=False) -> dict:
        """
        Get the summary for all of the samples

        If the summary has already been calculated, it will return the cached summary

        Args:
            reindex (bool): If True, reindex the summary DataFrames

        Returns:
            dict: A dictionary of DataFrames, where keys are parameter names and each
            DataFrame contains the summary statistics for a single parameter
        """
        if self.summary and not reindex:
            return self.summary

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
        stat_functions = {
            "median": np.median,
        }
        pyjags_data = az.from_pyjags(posterior=samples)
        summary = az.summary(pyjags_data, stat_funcs=stat_functions)
        for param_name, single_param_summary in summary.iterrows():
            if param_name not in samples_dict:
                samples_dict[param_name] = pd.DataFrame()
            single_param_summary_df = single_param_summary.to_frame().T
            samples_dict[param_name] = pd.concat(
                [samples_dict[param_name], single_param_summary_df], ignore_index=True
            )

    def save(self, directory: Union[Path, str]) -> None:
        """
        Save the results to a directory

        The summary of the results are stored in csv files with their name corresponding to
        the parameter.

        The params are stored in `metadata.json`.

        The raw samples are stored in `raw_samples/sim_n.npz` where each simulation
        gets its own file.

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

        summary_dict = self.get_summary(reindex=True)

        for param_name, summary_df in summary_dict.items():
            summary_df.to_csv(directory / f"{param_name}_summary.csv")

        raw_samples_dir = directory / "raw_samples"
        raw_samples_dir.mkdir(parents=True, exist_ok=True)

        for i, raw_sample in enumerate(self.samples_list):
            np.savez_compressed(
                file=raw_samples_dir / f"sim_{i}.npz",
                **raw_sample,
            )

    @staticmethod
    def load(directory: Union[str, Path]) -> "SimResults":
        """
        Load SimResults from directory

        You can save SimResults using the `save_summary` method

        Args:
            directory (str | Path)
        """
        directory = Path(directory)
        with open(directory / "metadata.json") as f:
            sim_params = json.load(f)

        # Get the allowed keys from the __init__ method of SimParams
        allowed_keys = inspect.signature(SimParams.__init__).parameters.keys()
        # Filter sim_params to only include allowed keys
        # some parameters are derived, and we can't pass those to the SimParams class
        filtered_sim_params = {k: v for k, v in sim_params.items() if k in allowed_keys}

        # get raw samples
        raw_samples_files = os.listdir(directory / "raw_samples")
        samples_list = []
        for raw_samples_file in raw_samples_files:
            raw_samples_file_path = directory / "raw_samples" / raw_samples_file
            sample = np.load(raw_samples_file_path)
            samples_list.append(sample)

        sim_results = SimResults(SimParams(**filtered_sim_params))
        sim_results.samples_list = samples_list

        # get the summary
        summary_files = os.listdir(directory)
        summary_files = [
            directory / summary_file for summary_file in summary_files if "summary" in summary_file
        ]
        summary_dict = {}
        for summary_file in summary_files:
            param_name = summary_file.stem.replace("_summary", "")
            summary_dict[param_name] = pd.read_csv(summary_file)
        sim_results.summary = summary_dict

        return sim_results

    def plot_summary_hists(self, params: Union[List[str], None] = None, func="median") -> None:
        """
        Makes a histogram for each parameter of the summary statistic specified by `func`

        Passing None to params will make a histogram for every parameter

        Args:
            params (List[str]): list of parameters to make histograms for
            func (str): the function to plot the histogram for. Default is "median"
        """
        if not self.summary:
            self.get_summary()

        summary = self.summary

        # Filter the parameters to be plotted
        derived = {"mean_psi", "PropOcc"}
        if params is None:
            params = [param for param in summary.keys() if param not in derived]
        else:
            params = [param for param in params if param in summary and param not in derived]
        
        print(summary.keys())

        num_params = len(params)
        num_cols = 3  # You can adjust the number of columns as needed
        num_rows = math.ceil(num_params / num_cols)

        _, axes = plt.subplots(num_rows, num_cols, figsize=(15, num_rows * 5))
        axes = axes.flatten()  # Flatten the axes array for easy iteration

        for i, param_name in enumerate(params):
            summary_df = summary[param_name]
            try:
                sim_param_value = self._get_sim_param_value(param_name)
            except ValueError:
                continue
            summary_values = summary_df[func]
            ax = axes[i]
            ax.axvline(x=sim_param_value, color="r", linestyle="--")
            ax.hist(summary_values, bins=20)
            ax.set_title(param_name)

        plt.tight_layout()

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
        print(sim_params, sim_param_name)
        if "[" not in sim_param_name and "]" not in sim_param_name and sim_param_name[-1] != "]":
            if sim_param_name not in sim_params:
                raise ValueError(f"param {sim_param_name} not found in sim_params")
            return sim_params[sim_param_name]
        name = "".join(sim_param_name.split("[")[0:-1])
        index = int(sim_param_name.split("[")[-1].replace("]", ""))
        try:
            value = sim_params[name][index]
        except ValueError:
            raise ValueError(f"name {name} and index {index} not found")
        return value
