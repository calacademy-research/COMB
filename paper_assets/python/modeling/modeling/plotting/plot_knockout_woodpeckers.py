from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from model_runs.model_runs_db import ResultsDB

db = ResultsDB.create(str(Path(__file__).parent.parent / "data" / "results_db"))

species = [
    "whhwoo",
    "rebsap",
    "pilwoo",
    "norfli",
    "haiwoo",
    "bkbwoo",
]

# species = ["lazbun", "macwar", "yerwar", "wewpew", "btywar", "herwar"]

models = {
    "single_year_single_species_all": "full",
    "single_year_single_species_no_pc": "nopc",
    "single_year_single_species_no_aru": "noaru",
    "single_year_single_species_no_scores": "noscores",
    "single_year_single_species_no_aru_no_pc": "scoresonly",
    "single_year_single_species_no_scores_no_pc": "aruonly",
    "single_year_single_species_no_scores_no_aru": "pconly",
}

# Collect data for all species and models
data = {}
for s in species:
    data[s] = {}
    for model_name, short_name in models.items():
        run_id = db.get_all_run_ids(dict(species=[s], aru_threshold=-0.5), model_name=model_name)
        print(s, model_name, run_id)
        run = db.get_run(run_id[0])
        mean_psi = run.results["posterior"].data_vars["mean_psi"].to_numpy().flatten()
        data[s][short_name] = mean_psi

# Create subplots: 2 rows, 3 columns
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

# Plot each species in its own subplot
for idx, s in enumerate(species):
    ax = axes[idx]

    # Prepare data for violin plot
    violin_data = [data[s][short_name] for short_name in models.values()]
    labels = list(models.values())
    positions = np.arange(1, len(labels) + 1)

    # Create violin plot
    vp = ax.violinplot(
        violin_data, positions=positions, widths=0.7, showmeans=True, showmedians=True
    )

    # Customize violin plot appearance
    for pc in vp["bodies"]:
        pc.set_facecolor("lightblue")
        pc.set_alpha(0.7)
        pc.set_edgecolor("black")
        pc.set_linewidth(1)

    ax.set_title(s.upper(), fontsize=12, fontweight="bold")
    ax.set_ylabel("Mean psi", fontsize=10)
    ax.set_xlabel("Model", fontsize=10)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")

plt.tight_layout()
plt.savefig("data/plots/knockout_woodpeckers_boxplot.pdf")
