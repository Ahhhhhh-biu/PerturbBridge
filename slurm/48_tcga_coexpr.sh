#!/bin/bash
#SBATCH --job-name=pb_coexpr
#SBATCH --partition=gpu_a100_80gb,gpu_a100_40gb,gpu_rtx8000_48gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/tcga_coexpr.txt
# RBP expression x event PSI co-variation across TCGA tumours (observational prior for the programme-level model).
PB=/well/blagden/users/eic789/PerturbBridge
$PB/envs/pb_gpu_venv/bin/python -W ignore -u $PB/13_model/programme/02_tcga_coexpr.py && echo "DONE coexpr"
