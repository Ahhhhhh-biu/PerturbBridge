#!/bin/bash
#SBATCH --job-name=pb_h4b2
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/h4_batch2_job_%x.log
# Inverse mapping on one cohort of the second batch (14_inverse/H4_BATCH2_PLAN.md). Usage: sbatch --job-name=pb_h4b2_<GSE> --export=ALL,GSE=<GSE> 53_h4_batch2.sh
PB=/well/blagden/users/eic789/PerturbBridge
[ -f $PB/05_precancer_cohorts/$GSE/rmats_all/.done ] || { echo "rMATS post not complete for $GSE"; exit 1; }
$PB/envs/pb_gpu_venv/bin/python -W ignore -u $PB/14_inverse/16_h4_batch2.py $GSE > $PB/logs/h4_batch2_$GSE.txt 2>&1 && echo "DONE h4b2 $GSE"
