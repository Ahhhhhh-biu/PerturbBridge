#!/bin/bash
#SBATCH --job-name=pb_tma
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/tma_job.log
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp
$PB/envs/pb_gpu_venv/bin/python $PB/14_inverse/08_tma_consistency.py > $PB/logs/tma_consistency.txt 2>&1 && echo "DONE tma $(date)"
