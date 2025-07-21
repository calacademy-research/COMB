import math
import os
import arviz as az
import pandas as pd
import polars as pl
import numpy as np
from .sim_data import SimParams, DataParams, ModelParams
from typing import List, Union
from pathlib import Path
import matplotlib.pyplot as plt
import json
import inspect
from tqdm import tqdm
from matplotlib.lines import Line2D


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

        self.within_ci = self.all_within_ci()

    def all_within_ci(self) -> bool:
        """
        Check if the true parameter values are within the 95% CI from the many samples of the posterior

        Returns:
            bool: True if all the true parameter values are within the 95% CI
        """
        if not self.summary:
            self.get_summary()

        for param_name, summary_df in self.summary.items():
            try:
                sim_param_value = self._get_sim_param_value(param_name)
            except ValueError:
                continue
            m = summary_df["mean"].mean()
            sd = summary_df["mean"].std()
            ci = 1.96 * sd
            lower = m - ci
            upper = m + ci
            if not (lower <= sim_param_value <= upper):
                return False
        return True

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
        if len(self.samples_list) == 0:
            raise ValueError(
                "No samples to save. Make sure that you have loaded the samples."
            )

        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)

        # make params into dict, and turn np arrays into lists for easier serialization
        metadata = {}
        metadata["model_params"] = self.sim_params.model_params.__dict__
        metadata["data_params"] = self.sim_params.data_params.__dict__
        for key, value in metadata["model_params"].items():
            if isinstance(value, np.ndarray):
                metadata["model_params"][key] = value.tolist()
        for key, value in metadata["data_params"].items():
            if isinstance(value, np.ndarray):
                metadata["data_params"][key] = value.tolist()

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
    def load(
        directory: Union[str, Path], load_raw_samples: bool = True
    ) -> "SimResults":
        """
        Load SimResults from directory

        You can save SimResults using the `save_summary` method

        Args:
            directory (str | Path)
            load_raw_samples (bool): If True, will load the raw samples. Default is True
        """
        directory = Path(directory)
        with open(directory / "metadata.json") as f:
            sim_params = json.load(f)

        # Get the allowed keys from the __init__ method of ModelParams and DataParams
        data_allowed_keys = inspect.signature(DataParams.__init__).parameters.keys()
        model_allowed_keys = inspect.signature(ModelParams.__init__).parameters.keys()
        # Filter sim_params to only include allowed keys
        # some parameters are derived, and we can't pass those to the SimParams class
        data_filtered_params = {
            k: v for k, v in sim_params.items() if k in data_allowed_keys
        }
        model_filtered_params = {
            k: v for k, v in sim_params.items() if k in model_allowed_keys
        }

        samples_list = []
        if load_raw_samples:
            # get raw samples
            raw_samples_files = os.listdir(directory / "raw_samples")
            for raw_samples_file in raw_samples_files:
                raw_samples_file_path = directory / "raw_samples" / raw_samples_file
                sample = np.load(raw_samples_file_path)
                samples_list.append(sample)

        sim_results = SimResults(
            SimParams(
                ModelParams(**model_filtered_params), DataParams(**data_filtered_params)
            )
        )
        sim_results.samples_list = samples_list

        # get the summary
        summary_files = os.listdir(directory)
        summary_files = [
            directory / summary_file
            for summary_file in summary_files
            if "summary" in summary_file
        ]
        summary_dict = {}
        for summary_file in summary_files:
            param_name = summary_file.stem.replace("_summary", "")
            summary_dict[param_name] = pd.read_csv(summary_file)
        sim_results.summary = summary_dict

        return sim_results

    def plot_summary_hists(
        self, params: Union[List[str], None] = None, func="median"
    ) -> None:
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
            params = [
                param for param in params if param in summary and param not in derived
            ]

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

            # compute the 95% credible interval (SDOM)
            sd = summary_df["mean"].std()
            m = summary_df["mean"].mean()
            ci = 1.96 * sd
            lower = m - ci
            upper = m + ci

            ax.errorbar(
                x=sim_param_value,
                y=max(np.histogram(summary_values, bins=20)[0]),
                xerr=[[m - lower], [upper - m]],
                fmt="o",
                color="r",
                capsize=5,
            )

            ax.axvline(x=sim_param_value, color="r", linestyle="--")
            ax.axvline(x=m, color="b", linestyle="--")
            ax.hist(summary_values, bins=20)
            ax.set_title(param_name)
            handles = [
                Line2D([0], [0], color="r", linestyle="--", label="Simulated Value"),
                Line2D([0], [0], color="b", linestyle="--", label="Mean Value"),
                Line2D([0], [0], color="r", marker="o", linestyle="", label="95% CI"),
            ]
            plt.legend(handles=handles, loc="upper right")

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
        if (
            "[" not in sim_param_name
            and "]" not in sim_param_name
            and sim_param_name[-1] != "]"
        ):
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


class ManySimResults:
    def __init__(self, sims_dir: Union[str, Path]):
        """
        Initialize a ManySimResults object

        Args:
            sims_dir (str | Path): The directory containing the simulation results
        """
        self.sims_dir = Path(sims_dir)
        self.sim_results: dict[str, SimResults] = {}
        self.load_sim_results()
        self.sim_results_index = self._create_sim_results_index()

    def load_sim_results(self):
        """
        Load all of the simulation results from the directory
        """
        sim_dirs = [x for x in self.sims_dir.iterdir() if x.is_dir()]
        for sim_dir in tqdm(sim_dirs):
            self.sim_results[sim_dir.name] = SimResults.load(
                sim_dir, load_raw_samples=False
            )

    def _create_sim_results_index(self):
        """
        Create an index for the simulation results

        Returns:
            pl.DataFrame: A DataFrame containing the metadata for each simulation
        """
        metadata = []
        for dir_name, sim_result in self.sim_results.items():
            params = sim_result.sim_params.__dict__
            params["dir_name"] = dir_name
            first_key = list(sim_result.summary.keys())[0]
            params["n_sims"] = len(sim_result.summary[first_key])
            params["within_ci"] = sim_result.within_ci
            metadata.append(params)
        return pl.DataFrame(metadata, strict=False)

    def load_from_sim_results_index(self, index: pl.DataFrame) -> List[SimResults]:
        """
        Load simulation results from the index.

        Ideally, this index is a subset of the `sim_results_index` attribute
        that contains simulations that are of interest.

        Args:
            index (pl.DataFrame): A DataFrame containing the metadata for each simulation

        Returns:
            List[SimResults]: A list of SimResults objects
        """
        sim_results: List[SimResults] = []

        for sim in index.iter_rows(named=True):
            dir_name = sim["dir_name"]
            sim_result = SimResults.load(
                self.sims_dir / dir_name, load_raw_samples=True
            )
            sim_results.append(sim_result)

        return sim_results
