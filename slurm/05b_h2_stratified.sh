#!/bin/bash
#SBATCH --job-name=pb_h2b
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/h2_stratified.log
PB=/well/blagden/users/eic789/PerturbBridge
$PB/envs/pb_venv/bin/python $PB/04_h2_check/02_h2_stratified.py
