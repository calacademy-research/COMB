import pandas as pd
import pyjags
import numpy as np
import pyjags.model
from sim_results import SimResults
from sim_data import SimParams
from sim_model import SimModel
import arviz as az

N_SIMS = 1

params = SimParams(
    nsites=25,
    nsurveys_aru=3,
    nsurveys_pc=3,
    covar_continuous=True,
    beta0=0,
    beta1=1,
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

summary = sim_results.get_summary()
print(summary)
