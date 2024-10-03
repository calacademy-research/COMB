import numpy as np
from typing import List, Tuple, Union
from pathlib import Path
from sim_lib.sim_data import SimParams
from sim_lib.sim_results import SimResults, ManySimResults
import os
from tqdm import tqdm

root_dir = Path("outputs/testing")
new_dir = Path("outputs/testing2")

dirs = os.listdir(root_dir)

for dir in tqdm(dirs):
    if not dir.startswith("sim_summary"):
        continue
    sim = SimResults.load(root_dir / dir)
    h = hash(sim.sim_params)
    existing_hashes = {int(n.split("_")[-1]) for n in os.listdir(new_dir)}
    if h not in existing_hashes:
        sim.save(new_dir / f"sim_summary_{h}")
    else:
        print(f"Loaded {h} because it already exists")
        loaded = SimResults.load(new_dir / f"sim_summary_{h}")
        for sample in sim.samples_list:
            loaded.append_samples(sample)
        loaded.save(new_dir / f"sim_summary_{h}")
