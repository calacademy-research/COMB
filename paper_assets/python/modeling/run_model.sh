#!/bin/bash

#SBATCH --job-name=run_caples_model
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --nodelist=flor
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err


uv run combined_model.py