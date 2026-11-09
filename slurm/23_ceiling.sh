#!/bin/bash
#SBATCH --job-name=pb_ceiling
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/ceiling_%j.log
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp
$PB/envs/pb_venv/bin/python $PB/12_baselines/01_ceiling.py && echo "DONE ceiling $(date)"
