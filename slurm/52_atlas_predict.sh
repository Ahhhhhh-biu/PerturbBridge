#!/bin/bash
#SBATCH --job-name=pb_atlas
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/atlas_predict.txt
# Predicted programme atlas for all human RBPs + HeyA8 LARP1 knockout check (13_model/programme/PLAN.md addendum b).
PB=/well/blagden/users/eic789/PerturbBridge
$PB/envs/pb_gpu_venv/bin/python -W ignore -u $PB/13_model/programme/09_atlas_predict.py && echo "DONE atlas"
