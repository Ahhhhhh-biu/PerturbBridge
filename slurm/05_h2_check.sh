#!/bin/bash
#SBATCH --job-name=pb_h2
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/h2_check.log
PB=/well/blagden/users/eic789/PerturbBridge
$PB/envs/pb_venv/bin/python $PB/04_h2_check/01_h2_perturbation_specific_variance.py
