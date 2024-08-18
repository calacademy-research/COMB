import pyjags
import numpy as np
import pyjags.model
from sim_data import SimParams, SimData
from sim_model import SimModelParams, SimModel

params = SimParams(
    nsites=10,
    nsurveys_aru=5,
    nsurveys_pc=5,
    covar_continuous=True,
).generate_sim_params()

sim_data = SimData(params)
sim_data.gen_simulated_data()

sim_model_params = SimModelParams(
    include_aru_model=True,
    include_scores_model=True,    
    sim_data=sim_data)
model = SimModel(sim_model_params)
samples = model.sample()
print(samples)