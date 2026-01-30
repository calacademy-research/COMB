#!/bin/bash

#SBATCH --job-name=run_comb_sim_model
#SBATCH --cpus-per-task=60
#SBATCH --mem=64G
#SBATCH --nodelist=rosalindf
#SBATCH --output=logs/%x-%j.out


uv run run_sim.py