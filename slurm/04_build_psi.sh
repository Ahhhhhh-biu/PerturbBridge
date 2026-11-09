#!/bin/bash
#SBATCH --job-name=pb_psi
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/build_psi.log
PB=/well/blagden/users/eic789/PerturbBridge
$PB/envs/pb_venv/bin/python $PB/03_suppa_psi/01_build_matrices_and_psi.py
