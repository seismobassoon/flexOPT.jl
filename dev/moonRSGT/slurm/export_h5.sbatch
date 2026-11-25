#!/bin/bash
# dev/moonRSGT/slurm/export_h5.sbatch -- JLD2 catalogue -> plain HDF5 (for Python / the students' hard disk)
#   cd /gpfs/users/nobuaki/github/flexOPT.jl && sbatch dev/moonRSGT/slurm/export_h5.sbatch
# Needs ≈ 0.6 TB more on scratch (the JLD2 catalogue can be deleted afterwards if only Python is used).
#SBATCH --job-name=moonRSGT-h5
#SBATCH --partition=ncpushort
#SBATCH --array=1-100%20
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/scratch/nobuaki/moonRSGT/logs/%x_%A_%a.out

export JULIA_DEPOT_PATH=/gpfs/users/nobuaki/.julia
export MOON_RSGT_DIR=/gpfs/scratch/nobuaki/moonRSGT/catalogue
export MOON_RSGT_H5_DIR=/gpfs/scratch/nobuaki/moonRSGT/catalogue_h5
conda deactivate 2>/dev/null || true
cd /gpfs/users/nobuaki/github/flexOPT.jl
julia --project=. dev/moonRSGT/export_h5.jl --task=${SLURM_ARRAY_TASK_ID} --per-task=10
