#!/bin/bash
#SBATCH --job-name=pb_h3tcga_b
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/h3_tcga_job_b.log
# Block 2b. First checks that all 8 Snaptron shards finished (logs/snaptron_<i>.log has "DONE <i> "); shards that timed out are resubmitted
# in resume mode (short partition, 6 h) and this job requeues itself behind them. Then runs 14_inverse/07_h3_tcga.py.
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; export OMP_NUM_THREADS=8
missing=""; for i in 0 1 2 3 4 5 6 7; do grep -q "^DONE $i " $PB/logs/snaptron_$i.log 2>/dev/null || missing="$missing,$i"; done; missing=${missing#,}
if [ -n "$missing" ]; then
  R=$(sbatch --parsable --export=ALL,NTASK=8 --array=$missing --partition=short --time=06:00:00 $PB/slurm/35_snaptron_tcga.sh)
  echo "$(date) shards $missing unfinished -> resumed as $R; requeue h3_tcga after it"
  N=$(sbatch --parsable --partition=short --time=06:00:00 --dependency=afterany:$R $PB/slurm/38_h3_tcga.sh); echo "$N pb_h3tcga requeued after resume array $R (shards $missing), $(date '+%Y-%m-%d %H:%M')" >> $PB/logs/JOBS.txt; echo $N > $PB/logs/.h3tcga_job; exit 0
fi
$PB/envs/pb_gpu_venv/bin/python -u $PB/14_inverse/07_h3_tcga.py > $PB/logs/h3_tcga_b.txt 2>&1 && echo "DONE h3_tcga $(date)"
