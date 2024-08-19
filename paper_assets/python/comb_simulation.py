import pandas as pd
import pyjags
import numpy as np
import pyjags.model
from sim_results import SimResults
from sim_data import SimParams
from sim_model import SimModel
import arviz as az

N_SIMS = 10

params = SimParams(
    nsites=80,
    nsurveys_aru=24,
    nsurveys_pc=3,
    covar_continuous=True,
    psi=0.1,
    threshold=0,
    include_covar_model=False,
)


def run_sim():
    model = SimModel(params)
    samples = model.sample()
    return samples


sim_results = SimResults(params)

for _ in range(N_SIMS):
    samples = run_sim()
    sim_results.append_samples(samples)

sim_results.save_summary("sim_summary")
