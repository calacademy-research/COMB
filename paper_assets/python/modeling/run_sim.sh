#!/bin/bash

#SBATCH --job-name=run_caples_sim
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --nodelist=flor
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err


uv run run_sim.py   