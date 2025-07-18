from typing import List, Tuple, Union
from pathlib import Path
from sim_lib.sim_data import SimParams
from sim_lib.sim_model import SimModel
from sim_lib.sim_results import SimResults
import itertools
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

N_SIMS = 99
N_PROCESSES = 16


class SimCombinations:
    def __init__(
        self,
        param_combinations: dict,
        output_dir: str = "sim_results",
        skip_existing: bool = False,
    ):
        """
        Initialize the SimCombinations object
        """
        self.skip_existing = skip_existing
        self.output_dir = Path(output_dir)
        self.param_combinations = param_combinations
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
        with ProcessPoolExecutor(max_workers=N_PROCESSES) as executor:
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

        for _ in range(N_SIMS):
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


if __name__ == "__main__":
    param_combinations = {
        "p11": [0.9, 0.5, 0.1],
        "p_aru11": [0.9, 0.5, 0.1],
        "p_aru01": [0.05, 0],
        "mu": [(-2, -1.75), (-2, 0)],
        "sigma": [(0.25, 1)],
        "threshold": [-1, 0, 1],
        "nsites": [80, 200],
        "include_covar_model": [True],
        "covar_continuous": [True],
        "beta0": [-1],
        "beta1": [-1, 1],
        "nsurveys_aru": [24],
        "nsurveys_scores": [24, 8],
        "nsurveys_pc": [3],
        "aru_scores_independent_model": [True, False],
    }
    combs = SimCombinations(
        param_combinations, output_dir="outputs/testing2", skip_existing=False
    )
    combs.run_all_combinations()
