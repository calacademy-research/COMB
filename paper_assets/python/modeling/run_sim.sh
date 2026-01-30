#!/bin/bash

#SBATCH --job-name=run_comb_sim
#SBATCH --cpus-per-task=62
#SBATCH --mem=64G
#SBATCH --nodelist=rosalindf
#SBATCH --output=logs/%x-%j.out


uv run run_sim.py   