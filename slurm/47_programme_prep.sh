#!/bin/bash
#SBATCH --job-name=pb_progprep
#SBATCH --partition=gpu_a100_80gb,gpu_a100_40gb,gpu_rtx8000_48gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/programme_prep.txt
# Matrices and RBP descriptor blocks for the programme-level model (13_model/programme/PLAN.md).
PB=/well/blagden/users/eic789/PerturbBridge; PY=$PB/envs/pb_gpu_venv/bin/python
$PY -W ignore -u $PB/13_model/programme/00_build_matrices.py && $PY -W ignore -u $PB/13_model/programme/01_rbp_descriptors.py && echo "DONE progprep"
