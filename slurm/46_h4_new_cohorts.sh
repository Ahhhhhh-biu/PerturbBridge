#!/bin/bash
#SBATCH --job-name=pb_h4new
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/h4_new_cohorts_job.log
# After the rMATS post jobs of GSE182338 and GSE289499: group-level summary, marker composition scores, then 14_inverse/14_h4_new_cohorts.py.
PB=/well/blagden/users/eic789/PerturbBridge; PY=$PB/envs/pb_gpu_venv/bin/python; export OMP_NUM_THREADS=8; export TMPDIR=$PB/envs/tmp
for pair in "GSE182338 rmats_carrier_vs_rest" "GSE289499 rmats_progression"; do set -- $pair; D=$PB/05_precancer_cohorts/$1/$2
  [ -f $D/.done ] || { echo "post not done for $1 $2"; exit 1; }
  $PY $PB/05_precancer_cohorts/summarise_rmats.py $D $1 > $D/summary.txt 2>&1 && $PY $PB/05_precancer_cohorts/cell_composition_scores.py $PB/05_precancer_cohorts/$1 $D > $PB/logs/composition_$1.txt 2>&1 || echo "composition step failed for $1 (see logs/composition_$1.txt)"
done
$PY -u $PB/14_inverse/14_h4_new_cohorts.py > $PB/logs/h4_new_cohorts.txt 2>&1 && echo "DONE h4new $(date)"
