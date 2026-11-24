#!/bin/bash
# dev/moonRSGT/slurm/submit_all.sh -- the whole 1000-model catalogue in one command (run on stella, it only submits).
#
#   cd /gpfs/users/nobuaki/github/flexOPT.jl && bash dev/moonRSGT/slurm/submit_all.sh
#
# 1. prepare + check (model 00001, rebuilt with the current parameters, + comparison figure)
# 2. the array (100 tasks × 10 models, ≤ 50 at a time), which starts ONLY if step 1 succeeded.
#    Look at the check log/figure while the array waits or runs; if it looks wrong: scancel <array jobid>.
# Re-running this script later is safe: finished models are skipped.
set -euo pipefail
cd "$(dirname "$0")/../../.."
mkdir -p /gpfs/scratch/nobuaki/moonRSGT/logs /gpfs/scratch/nobuaki/moonRSGT/catalogue
ls /gpfs/scratch/nobuaki/moonRSGT/OUTPUT_ND/*.nd > /dev/null    # the confidential models must be there
jid=$(sbatch --parsable dev/moonRSGT/slurm/prepare_check.sbatch)
aid=$(sbatch --parsable --dependency=afterok:${jid} dev/moonRSGT/slurm/build_catalogue.sbatch)
echo "check job  : ${jid}  -> /gpfs/scratch/nobuaki/moonRSGT/logs/moonRSGT-check_${jid}.out"
echo "            figure -> dev/moonRSGT/outputs/check/check_00001.png"
echo "array job  : ${aid}  (starts after ${jid} succeeds) -> logs/moonRSGT_${aid}_<task>.out"
echo "follow     : squeue -u \$USER ; seff ${jid}"
echo "afterwards : srun -p ncpushort -n1 -c1 --mem=4G -t 0:30:00 julia --project=. dev/moonRSGT/run.jl info"
