#!/bin/bash
# dev/moonRSGT/slurm/build_catalogue.sbatch -- Moon RSGT catalogue, one Slurm array task per 10 models.
#
#   mkdir -p /gpfs/scratch/nobuaki/moonRSGT/logs
#   cd /gpfs/users/nobuaki/github/flexOPT.jl
#   julia --project=. dev/moonRSGT/run.jl prepare          # once, BEFORE the array (caches the OPT3 recipes)
#   sbatch dev/moonRSGT/slurm/build_catalogue.sbatch
#   julia --project=. dev/moonRSGT/run.jl info             # afterwards: what is missing
#
# 1000 models / 10 per task = 100 tasks; at most 50 run at the same time (%50).
# Re-submitting is safe: finished models are skipped (same parameters).
# Adjust --time after timing one model with `run.jl check`.
#SBATCH --job-name=moonRSGT
#SBATCH --array=1-100%50
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/scratch/nobuaki/moonRSGT/logs/%x_%A_%a.out
##SBATCH --partition=CHANGE_ME          # remove one '#' and set stella's partition if needed

export JULIA_DEPOT_PATH=/gpfs/users/nobuaki/.julia
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MOON_ND_DIR=/gpfs/scratch/nobuaki/moonRSGT/OUTPUT_ND
export MOON_RSGT_DIR=/gpfs/scratch/nobuaki/moonRSGT/catalogue

cd /gpfs/users/nobuaki/github/flexOPT.jl
julia --project=. dev/moonRSGT/run.jl build --task=${SLURM_ARRAY_TASK_ID} --per-task=10
