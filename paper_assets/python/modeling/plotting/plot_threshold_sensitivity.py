from pathlib import Path
import matplotlib.pyplot as plt

import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))

from model_runs.model_runs_db import ResultsDB

db = ResultsDB.create(str(Path(__file__).parent.parent / "data" / "results_db"))

species = "herwar"

threshold_values = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]

# Collect data for all thresholds
p_aru11_data = []
p_aru01_data = []

for threshold in threshold_values:
    run_id = db.get_all_run_ids(dict(species=[species], aru_threshold=threshold))
    run = db.get_run(run_id[0])

    p_aru11_values = run.results["posterior"].data_vars["p_aru11"].to_numpy().flatten()
    p_aru01_values = run.results["posterior"].data_vars["p_aru01"].to_numpy().flatten()

    p_aru11_data.append(p_aru11_values)
    p_aru01_data.append(p_aru01_values)

# Create side-by-side box plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot p_aru11
bp1 = ax1.boxplot(
    p_aru11_data, tick_labels=[str(t) for t in threshold_values], patch_artist=True
)
ax1.set_xlabel("Threshold", fontsize=12)
ax1.set_ylabel("p_aru11", fontsize=12)
ax1.set_title("p_aru11 Distribution by Threshold", fontsize=14, fontweight="bold")
ax1.grid(True, alpha=0.3)

# Color the boxes
for patch in bp1["boxes"]:
    patch.set_facecolor("lightblue")
    patch.set_alpha(0.7)

# Plot p_aru01
bp2 = ax2.boxplot(
    p_aru01_data, tick_labels=[str(t) for t in threshold_values], patch_artist=True
)
ax2.set_xlabel("Threshold", fontsize=12)
ax2.set_ylabel("p_aru01", fontsize=12)
ax2.set_title("p_aru01 Distribution by Threshold", fontsize=14, fontweight="bold")
ax2.grid(True, alpha=0.3)

# Color the boxes
for patch in bp2["boxes"]:
    patch.set_facecolor("lightcoral")
    patch.set_alpha(0.7)

# Add overall title with species
fig.suptitle(
    f"Threshold Sensitivity Analysis - {species.upper()}",
    fontsize=16,
    fontweight="bold",
    y=0.98,
)

plt.tight_layout()
plt.savefig("data/plots/threshold_sensitivity_herwar.pdf")
