from typing import List, Tuple, Union
from pathlib import Path
from sim_lib.sim_data import ModelParams, SimParams, DataParams
from sim_lib.sim_model import SimModel, SimData
from sim_lib.sim_results import SimResults
import itertools
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
from sim_lib.sim_configs import configs
import os

DEFAULT_SIMS = 100
DEFAULT_PROCESSES = 16


class SimCombinations:
    def __init__(
        self,
        param_combinations: dict,
        output_dir: str = "sim_results",
        skip_existing: bool = False,
        n_sims: int = DEFAULT_SIMS,
        n_processes: int = DEFAULT_PROCESSES,
    ):
        """
        Initialize the SimCombinations object
        """
        # Split the model and data params dicts
        model_params_dict = param_combinations.get("model_params", {})
        data_params_dict = param_combinations.get("data_params", {})

        # some class vars
        self.skip_existing = skip_existing
        self.output_dir = Path(output_dir)
        self.n_sims = n_sims
        self.n_processes = n_processes

        # get the data params
        data_combinations = list(self.product_dict(*data_params_dict))
        self.valid_data_params = self.check_and_create_valid_data_params(
            data_combinations
        )
        print(
            f"Number of valid data parameter combinations: {len(self.valid_data_params)}"
        )

        # get the model params
        model_combinations = list(self.product_dict(*model_params_dict))
        self.valid_model_params = self.check_and_creat_valid_model_params(
            model_combinations
        )
        print(
            f"Number of valid model parameter combinations: {len(self.valid_data_params)}"
        )

    def product_dict(self, **kwargs):
        """
        Create a generator that yields all possible combinations of the input parameters

        Essentially computes the Cartesian product of the input parameters
        """
        keys = kwargs.keys()
        for instance in itertools.product(*kwargs.values()):
            yield dict(zip(keys, instance))

    def check_and_creat_valid_model_params(self, model_param_combinations: List[dict]):
        combinations: List[ModelParams] = []
        for comb in tqdm(model_param_combinations):
            model_params = ModelParams(**comb)
            combinations.append(model_params)

        return combinations

    def check_and_create_valid_data_params(self, data_param_combinations: List[dict]):
        """
        Given a list of data parameter combinations, check if each combination is valid
        """
        combinations: List[DataParams] = []
        for comb in tqdm(data_param_combinations):
            sim_params, is_valid = self._check_single_combination(comb)
            if is_valid and sim_params:
                combinations.append(sim_params)
        return combinations

    def _check_single_combination(
        self, combination: dict
    ) -> Tuple[Union[DataParams, None], bool]:
        """
        Check if a single combination of data parameters is valid

        Returns:
            Tuple[Union[DataParams, None], bool]: A tuple containing the DataParams object if the combination is valid,
            and a boolean indicating whether the combination is valid
        """
        # First uses the SimParams built-in __init__ method to check if the combination is valid
        # this will not catch every case, but it is a good starting point
        try:
            data_params = DataParams(**combination)
            model_params = ModelParams()  # Use default model params
            params = SimParams(data_params=data_params, model_params=model_params)
        except Exception:
            return None, False

        if not params.model_params.aru_scores_independent_model:
            if params.data_params.nsurveys_aru != params.data_params.nsurveys_scores:
                return None, False

        # we return the data params, as the model params are separate
        return data_params, True

    def run_all_combinations(self, dry_run: bool = False):
        """
        Run all of the valid parameter combinations using ProcessPoolExecutor

        Args:
            dry_run (bool): If True, will not run the simulations, but will print when they would have been run
        """
        with ProcessPoolExecutor(max_workers=self.n_processes) as executor:
            futures = [
                executor.submit(self.run_single_combination, params, dry_run)
                for params in self.valid_data_params
            ]
            for future in tqdm(as_completed(futures), total=len(futures)):
                try:
                    future.result()  # This will raise an exception if the task failed
                except Exception as e:
                    print(f"Simulation failed with exception: {e}")

    def run_single_combination(self, data_params: DataParams, dry_run: bool = False):
        """
        Run a single simulation given the parameters
        """
        # We need to run the data params across all models.
        # First we need to load all existing sim results into a map
        # where the key is the hash

        existing_dirs = [x for x in self.output_dir.iterdir() if x.is_dir()]
        existing_hashes = {int(x.name.split("_")[-1]) for x in existing_dirs}

        results_map: dict[int, SimResults] = {}

        for model_params in self.valid_model_params:
            params = SimParams(model_params, data_params)
            sim_params_hash = hash(params)
            if sim_params_hash in existing_hashes:
                if self.skip_existing:
                    print(f"Skipping existing simulation with parameters: {params}")
                    return
                sim_results = SimResults.load(
                    self.output_dir / f"sim_summary_{sim_params_hash}"
                )
            else:
                sim_results = SimResults(params)

            if dry_run:
                print(f"Would have run simulation with parameters: {params}")
                return

            results_map[sim_params_hash] = sim_results

        # TODO: for each sim, run across all models with same data

        for _ in range(self.n_sims):
            # first we generate data
            try:
                data_model = SimData(data_params)
                raw_sim_data, covars = data_model.get_data()
            except Exception as e:
                print(f"Error in generating data: {e}")
                print(f"Parameters: {data_params}")
                return

            # we fit the same data to each model
            for h, res in results_map.items():
                model_params = res.sim_params.model_params
                try:
                    model = SimModel(model_params, raw_sim_data, covars)
                    samples = model.sample()
                except Exception as e:
                    print(f"Error in running simulation: {e}")
                    print(f"Model parameters: {model_params}")
                    return

                results_map[h].append_samples(samples)

        # after we go through all models and all sims, write to disk
        for h, res in results_map.items():
            filepath = self.output_dir / f"sim_summary_{h}"
            res.save(filepath)


def parse_args():
    parser = argparse.ArgumentParser(description="Run simulation with specified config")
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Name of the configuration to use from sim_configs.py",
    )
    parser.add_argument(
        "--num_sims",
        type=int,
        default=DEFAULT_SIMS,
        help=f"Number of simulations to run (default: {DEFAULT_SIMS})",
    )
    parser.add_argument(
        "--num_processes",
        type=int,
        default=DEFAULT_PROCESSES,
        help=f"Number of processes to use (default: {DEFAULT_PROCESSES})",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    cfg = configs.get(args.config)
    if cfg is None:
        raise ValueError(f"config name {args.config} does not exist in sim_configs!")

    # name of output dir is the config name
    os.makedirs(f"data/{args.config}", exist_ok=True)

    combs = SimCombinations(
        cfg,
        output_dir=f"data/{args.config}",
        skip_existing=False,
        n_sims=args.num_sims,
        n_processes=args.num_processes,
    )
    combs.run_all_combinations()
