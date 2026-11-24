#!/bin/bash
# dev/moonRSGT/slurm/build_catalogue.sbatch -- Moon RSGT catalogue, one Slurm array task per 10 models.
#
#   mkdir -p /gpfs/scratch/nobuaki/moonRSGT/logs
#   cd /gpfs/users/nobuaki/github/flexOPT.jl
#   sbatch dev/moonRSGT/slurm/prepare_check.sbatch          # once, BEFORE the array (recipes cached + check)
#   sbatch dev/moonRSGT/slurm/build_catalogue.sbatch
#   srun -p ncpushort -n1 -c1 --mem=4G -t 0:30:00 julia --project=. dev/moonRSGT/run.jl info   # what is missing
#
# S-CAPAD/DANTE rules: never run julia on stella* (login), always give --mem
# (default 1929 MB/core on ncpu, the job is killed above it), keep %A/%a/%j in log names.
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
#SBATCH --partition=ncpu                 # 2 GB/core nodes, MaxTime 2 days

export JULIA_DEPOT_PATH=/gpfs/users/nobuaki/.julia
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MOON_ND_DIR=/gpfs/scratch/nobuaki/moonRSGT/OUTPUT_ND
export MOON_RSGT_DIR=/gpfs/scratch/nobuaki/moonRSGT/catalogue
conda deactivate 2>/dev/null || true

cd /gpfs/users/nobuaki/github/flexOPT.jl
julia --project=. dev/moonRSGT/run.jl build --task=${SLURM_ARRAY_TASK_ID} --per-task=10
