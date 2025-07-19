from typing import List, Tuple, Union
from pathlib import Path
from sim_lib.sim_data import SimParams
from sim_lib.sim_model import SimModel
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
        self.skip_existing = skip_existing
        self.output_dir = Path(output_dir)
        self.param_combinations = param_combinations
        self.n_sims = n_sims
        self.n_processes = n_processes
        self.combinations = list(self.product_dict(**param_combinations))
        self.valid_params = self.check_and_create_valid_params(self.combinations)
        print(f"Number of valid parameter combinations: {len(self.valid_params)}")

    def product_dict(self, **kwargs):
        """
        Create a generator that yields all possible combinations of the input parameters

        Essentially computes the Cartesian product of the input parameters
        """
        keys = kwargs.keys()
        for instance in itertools.product(*kwargs.values()):
            yield dict(zip(keys, instance))

    def check_and_create_valid_params(self, param_combinations: List[dict]):
        """
        Given a list of parameter combinations, check if each combination is valid
        """
        combinations: List[SimParams] = []
        for comb in tqdm(param_combinations):
            sim_params, is_valid = self._check_single_combination(comb)
            if is_valid and sim_params:
                combinations.append(sim_params)
        return combinations

    def _check_single_combination(
        self, combination: dict
    ) -> Tuple[Union[SimParams, None], bool]:
        """
        Check if a single combination of parameters is valid

        Returns:
            Tuple[Union[SimParams, None], bool]: A tuple containing the SimParams object if the combination is valid,
            and a boolean indicating whether the combination is valid
        """
        # First uses the SimParams built-in __init__ method to check if the combination is valid
        # this will not catch every case, but it is a good starting point
        try:
            params = SimParams(**combination)
        except Exception:
            return None, False

        if not params.aru_scores_independent_model:
            if params.nsurveys_aru != params.nsurveys_scores:
                return None, False

        return params, True

    def run_all_combinations(self, dry_run: bool = False):
        """
        Run all of the valid parameter combinations using ProcessPoolExecutor

        Args:
            dry_run (bool): If True, will not run the simulations, but will print when they would have been run
        """
        with ProcessPoolExecutor(max_workers=self.n_processes) as executor:
            futures = [
                executor.submit(self.run_single_combination, params, dry_run)
                for params in self.valid_params
            ]
            for future in tqdm(as_completed(futures), total=len(futures)):
                try:
                    future.result()  # This will raise an exception if the task failed
                except Exception as e:
                    print(f"Simulation failed with exception: {e}")

    def run_single_combination(self, params: SimParams, dry_run: bool = False):
        """
        Run a single simulation given the parameters

        Args:
            params (SimParams): The parameters to run the simulation with
            dry_run (bool): If True, will not run the simulation, but will print when it would have been run
        """
        existing_dirs = [x for x in self.output_dir.iterdir() if x.is_dir()]
        existing_hashes = {int(x.name.split("_")[-1]) for x in existing_dirs}
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

        for _ in range(self.n_sims):
            try:
                model = SimModel(params)
                samples = model.sample()
            except Exception as e:
                print(f"Error in running simulation: {e}")
                print(f"Parameters: {params}")
                return
            sim_results.append_samples(samples)

        filepath = self.output_dir / f"sim_summary_{sim_params_hash}"
        sim_results.save(filepath)


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
