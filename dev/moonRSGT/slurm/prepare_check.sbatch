#!/bin/bash
# dev/moonRSGT/slurm/prepare_check.sbatch -- `prepare` then `check` on a compute node
# (S-CAPAD/DANTE rule: nothing runs on the stella* login nodes).
#
#   mkdir -p /gpfs/scratch/nobuaki/moonRSGT/logs
#   cd /gpfs/users/nobuaki/github/flexOPT.jl
#   sbatch dev/moonRSGT/slurm/prepare_check.sbatch
#   squeue -u $USER ;  less /gpfs/scratch/nobuaki/moonRSGT/logs/moonRSGT-check_<jobid>.out
#   seff <jobid>      # memory / time actually used -> adjust build_catalogue.sbatch
#SBATCH --job-name=moonRSGT-check
#SBATCH --partition=ncpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/scratch/nobuaki/moonRSGT/logs/%x_%j.out

export JULIA_DEPOT_PATH=/gpfs/users/nobuaki/.julia
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MOON_ND_DIR=/gpfs/scratch/nobuaki/moonRSGT/OUTPUT_ND
export MOON_RSGT_DIR=/gpfs/scratch/nobuaki/moonRSGT/catalogue
conda deactivate 2>/dev/null || true

echo "Jobid $SLURM_JOB_ID on $SLURM_NODELIST"
cd /gpfs/users/nobuaki/github/flexOPT.jl
julia --project=. dev/moonRSGT/run.jl prepare && \
julia --project=. dev/moonRSGT/run.jl check
