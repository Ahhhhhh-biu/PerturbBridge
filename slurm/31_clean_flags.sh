#!/bin/bash
#SBATCH --job-name=pb_cleanflag
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/clean_flags.log
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/08_clean_flags_and_ceiling.py && echo "DONE clean_flags $(date)"
