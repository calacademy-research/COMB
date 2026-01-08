from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from model_runs.model_runs_db import ResultsDB

db = ResultsDB.create(str(Path(__file__).parent.parent / "data" / "results_db"))


species = ["herwar", "stejay", "yerwar", "wewpew", "mouqua", "btywar"]
model_name = "single_year_single_species_all"

# Create a 2x3 subplot layout
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for idx, s in enumerate(species):
    run_id = db.get_all_run_ids(dict(species=[s], aru_threshold=0), model_name)
    print(s, run_id)
    run = db.get_run(run_id[0])

    # Get posterior samples for mu and sigma
    mu0 = run.results["posterior"].data_vars["mu"].to_numpy().reshape(-1, 2)[:, 0]
    sigma0 = run.results["posterior"].data_vars["sigma"].to_numpy().reshape(-1, 2)[:, 0]

    mu1 = run.results["posterior"].data_vars["mu"].to_numpy().reshape(-1, 2)[:, 1]
    sigma1 = run.results["posterior"].data_vars["sigma"].to_numpy().reshape(-1, 2)[:, 1]

    # Use posterior means for plotting
    mu0_mean = np.mean(mu0)
    sigma0_mean = np.mean(sigma0)
    mu1_mean = np.mean(mu1)
    sigma1_mean = np.mean(sigma1)

    # Create x range for plotting
    x_min = min(mu0_mean - 3 * sigma0_mean, mu1_mean - 3 * sigma1_mean)
    x_max = max(mu0_mean + 3 * sigma0_mean, mu1_mean + 3 * sigma1_mean)
    x = np.linspace(x_min, x_max, 1000)

    # Plot distributions
    ax = axes[idx]
    ax.plot(
        x,
        norm.pdf(x, mu0_mean, sigma0_mean),
        label="Not Present",
        color="blue",
        linewidth=2,
    )
    ax.plot(
        x, norm.pdf(x, mu1_mean, sigma1_mean), label="Present", color="red", linewidth=2
    )

    ax.set_title(s.upper(), fontsize=12, fontweight="bold")
    ax.set_xlabel("Score", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("data/plots/score_distributions_by_species.pdf")
