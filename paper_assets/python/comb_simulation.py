import pyjags
import numpy as np
import pyjags.model
from sim_data import SimParams, SimData
from sim_model import SimModelParams, SimModel
import arviz as az

params = SimParams(
    nsites=25,
    nsurveys_aru=3,
    nsurveys_pc=3,
    covar_continuous=True,
    beta0=0,
    beta1=1,
)

sim_data = SimData(params)
sim_data.gen_simulated_data()

sim_model_params = SimModelParams(
    include_aru_model=True, include_scores_model=True, sim_data=sim_data
)
model = SimModel(sim_model_params)
samples = model.sample()

pyjags_data = az.from_pyjags(posterior=samples)
summary = az.summary(pyjags_data)
