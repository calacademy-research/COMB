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
import logging
import multiprocessing
from logging.handlers import QueueHandler, QueueListener
import atexit

DEFAULT_SIMS = 100
DEFAULT_PROCESSES = 16


def setup_logging(log_file: str):
    """
    Set up centralized logging for multiprocessing
    """
    # Create a queue for logging
    log_queue = multiprocessing.Queue()

    # Set up the listener process
    listener = QueueListener(
        log_queue,
        logging.FileHandler(log_file),
    )
    listener.start()

    # Configure the root logger to use the queue
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    root_logger.addHandler(QueueHandler(log_queue))

    # Set up formatter
    formatter = logging.Formatter(
        "%(asctime)s - PID:%(process)d - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Apply formatter to handlers
    for handler in listener.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.setFormatter(formatter)

    # Stop listener when program exits
    atexit.register(listener.stop)

    return log_queue, listener


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
        data_combinations = list(self.product_dict(**data_params_dict))
        self.valid_data_params = self.check_and_create_valid_data_params(
            data_combinations
        )
        logging.info(
            f"Number of valid data parameter combinations: {len(self.valid_data_params)}"
        )

        # get the model params
        model_combinations = list(self.product_dict(**model_params_dict))
        self.valid_model_params = self.check_and_create_valid_model_params(
            model_combinations
        )
        logging.info(
            f"Number of valid model parameter combinations: {len(self.valid_model_params)}"
        )

    def product_dict(self, **kwargs):
        """
        Create a generator that yields all possible combinations of the input parameters

        Essentially computes the Cartesian product of the input parameters
        """
        keys = kwargs.keys()
        for instance in itertools.product(*kwargs.values()):
            yield dict(zip(keys, instance))

    def check_and_create_valid_model_params(self, model_param_combinations: List[dict]):
        combinations: List[ModelParams] = []
        for comb in tqdm(model_param_combinations):
            model_params = ModelParams(**comb)
            if not (
                model_params.include_aru_model
                or model_params.include_pc_model
                or model_params.include_scores_model
            ):
                continue
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
        logging.info(
            f"Starting simulation with {len(self.valid_data_params)} data combinations"
        )
        logging.info(
            f"Using {self.n_processes} processes, {self.n_sims} simulations each"
        )

        with ProcessPoolExecutor(max_workers=self.n_processes) as executor:
            futures = [
                executor.submit(self.run_single_combination, params, dry_run)
                for params in self.valid_data_params
            ]
            completed = 0
            for future in tqdm(as_completed(futures), total=len(futures)):
                try:
                    future.result()  # This will raise an exception if the task failed
                    completed += 1
                    logging.info(
                        f"Completed {completed}/{len(futures)} parameter combinations"
                    )
                except Exception as e:
                    logging.error(f"Simulation failed with exception: {e}")

        logging.info("All simulations completed")

    def run_single_combination(self, data_params: DataParams, dry_run: bool = False):
        """
        Run a single simulation given the parameters
        """
        logging.info(f"Starting combination: {data_params}")

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
                    logging.info(
                        f"Skipping existing simulation with hash: {sim_params_hash}"
                    )
                    return
                sim_results = SimResults.load(
                    self.output_dir / f"sim_summary_{sim_params_hash}"
                )
                logging.info(f"Loaded existing results for hash: {sim_params_hash}")
            else:
                sim_results = SimResults(params)
                logging.info(f"Created new results for hash: {sim_params_hash}")

            if dry_run:
                logging.info(
                    f"DRY RUN: Would run simulation with hash: {sim_params_hash}"
                )
                return

            results_map[sim_params_hash] = sim_results

        for sim_num in range(self.n_sims):
            # first we generate data
            try:
                data_model = SimData(data_params)
                raw_sim_data, covars = data_model.get_data()
                logging.debug(
                    f"Generated data for simulation {sim_num + 1}/{self.n_sims}"
                )
            except Exception as e:
                logging.error(f"Error in generating data for sim {sim_num + 1}: {e}")
                logging.error(f"Parameters: {data_params}")
                return

            # we fit the same data to each model
            for h, res in results_map.items():
                try:
                    model = SimModel(res.sim_params, raw_sim_data, covars)
                    samples = model.sample()
                    results_map[h].append_samples(samples)
                    logging.debug(f"Completed model {h} for simulation {sim_num + 1}")
                except Exception as e:
                    logging.error(
                        f"Error in running simulation {sim_num + 1} for model {h}: {e}"
                    )
                    logging.error(f"Parameters: {res.sim_params}")
                    return

        # after we go through all models and all sims, write to disk
        for h, res in results_map.items():
            filepath = self.output_dir / f"sim_summary_{h}"
            res.save(filepath)
            logging.info(f"Saved results for model hash: {h}")

        logging.info(f"Completed all simulations for data combination: {data_params}")


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

    # Set up logging
    log_file = f"data/{args.config}_simulation.log"
    log_queue, listener = setup_logging(log_file)

    logging.info(f"Starting simulation run with config: {args.config}")
    logging.info(f"Log file: {log_file}")

    cfg = configs.get(args.config)
    if cfg is None:
        logging.error(f"config name {args.config} does not exist in sim_configs!")
        raise ValueError(f"config name {args.config} does not exist in sim_configs!")

    # name of output dir is the config name
    os.makedirs(f"data/{args.config}", exist_ok=True)
    logging.info(f"Created output directory: data/{args.config}")

    combs = SimCombinations(
        cfg,
        output_dir=f"data/{args.config}",
        skip_existing=False,
        n_sims=args.num_sims,
        n_processes=args.num_processes,
    )
    combs.run_all_combinations()

    logging.info("Simulation run completed")
