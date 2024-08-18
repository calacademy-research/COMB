import pandas as pd
import pyjags
import numpy as np
import pyjags.model
from sim_results import SimResults
from sim_data import SimParams, SimData
from sim_model import SimModelParams, SimModel
import arviz as az

N_SIMS = 2


def run_sim():
    params = SimParams(
        nsites=25,
        nsurveys_aru=3,
        nsurveys_pc=3,
        covar_continuous=True,
        beta0=0,
        beta1=1,
        threshold=-2,
    )
    sim_data = SimData(params)
    sim_data.gen_simulated_data()

    sim_model_params = SimModelParams(
        include_aru_model=True, include_scores_model=True, sim_data=sim_data
    )

    model = SimModel(sim_model_params)
    
    samples = model.sample()
    return samples


sim_results = SimResults()

for _ in range(N_SIMS):
    samples = run_sim()
    sim_results.append_samples(samples)

summary = sim_results.get_summary()
print(summary)
