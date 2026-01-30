from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from matplotlib.lines import Line2D

import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from model_runs.model_runs_db import ResultsDB

db = ResultsDB.create(str(Path(__file__).parent.parent / "data" / "results_db"))


species = ["herwar", "stejay", "yerwar", "wewpew", "mouqua", "btywar"]
model_name = "single_year_single_species_all_datetime"

# Collect alpha1 values for all species
alpha1_data = []
species_labels = []

for s in species:
    run_id = db.get_all_run_ids(dict(species=[s], aru_threshold=0), model_name)
    print(run_id, s)
    run = db.get_run(run_id[0])

    alpha1 = run.results["posterior"].data_vars["alpha1"].to_numpy().flatten()
    alpha1_data.append(alpha1)
    species_labels.append(s.upper())

# Create violin plot
fig, ax = plt.subplots(figsize=(12, 6))

parts = ax.violinplot(
    alpha1_data,
    positions=range(len(species)),
    showmeans=False,
    showmedians=True,
    widths=0.7,
)

# Customize the violin plot
for pc in parts["bodies"]:  # type: ignore
    pc.set_facecolor("#8dd3c7")
    pc.set_alpha(0.7)
    pc.set_edgecolor("black")
    pc.set_linewidth(1.5)

# Customize other elements
parts["cmedians"].set_color("blue")
parts["cmedians"].set_linewidth(2)

# Set labels and title
ax.set_xticks(range(len(species)))
ax.set_xticklabels(species_labels, rotation=45, ha="right")
ax.set_ylabel("Alpha1 (Time of Day Effect)", fontsize=12, fontweight="bold")
ax.set_xlabel("Species", fontsize=12, fontweight="bold")
ax.set_title(
    "Distribution of Alpha1 Parameter Across Species", fontsize=14, fontweight="bold"
)
ax.grid(axis="y", alpha=0.3, linestyle="--")

plt.tight_layout()
plt.savefig("data/plots/alpha1_violin_plot.pdf")
