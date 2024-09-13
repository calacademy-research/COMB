import pyjags.model
from sim_results import SimResults
from sim_data import SimParams
from sim_model import SimModel
from tqdm import tqdm

N_SIMS = 100

params = SimParams(
    nsites=80,
    nsurveys_aru=24,
    nsurveys_pc=3,
    covar_continuous=True,
    beta0=0.2,
    beta1=0.5,
    threshold=0,
    include_covar_model=True,
    sigma=(0.5, 2),
)


def run_sim():
    model = SimModel(params)
    samples = model.sample()
    return samples


sim_results = SimResults(params)

for _ in tqdm(range(N_SIMS)):
    samples = run_sim()
    sim_results.append_samples(samples)

sim_results.save("sim_summary")
