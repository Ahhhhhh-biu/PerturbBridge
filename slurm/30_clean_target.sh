#!/bin/bash
#SBATCH --job-name=pb_clean
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/clean_target.log
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/07_clean_target.py && echo "DONE clean_target $(date)"
