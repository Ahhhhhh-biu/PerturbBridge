#!/bin/bash
#SBATCH --job-name=pb_ceiling
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/ceiling.log
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/05_split_half_ceiling.py && echo "DONE ceiling $(date)"
