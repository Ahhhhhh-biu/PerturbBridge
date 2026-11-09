#!/bin/bash
#SBATCH --job-name=pb_h2c
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/h2_proper_ceiling.log
PB=/well/blagden/users/eic789/PerturbBridge
$PB/envs/pb_venv/bin/python $PB/04_h2_check/03_h2_proper_ceiling.py
