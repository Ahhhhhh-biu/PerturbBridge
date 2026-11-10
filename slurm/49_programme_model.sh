#!/bin/bash
#SBATCH --job-name=pb_prog
#SBATCH --partition=gpu_a100_80gb,gpu_a100_40gb,gpu_rtx8000_48gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/programme_model_%j.txt
# Programme-level forward model with permuted controls. Usage: sbatch [--export=ALL,TAG=_xxx,EXTRA=--no_cx] 49_programme_model.sh
PB=/well/blagden/users/eic789/PerturbBridge
$PB/envs/pb_gpu_venv/bin/python -W ignore -u $PB/13_model/programme/03_programme_model.py --tag "${TAG:-}" ${EXTRA:-} && $PB/envs/pb_gpu_venv/bin/python -W ignore -u $PB/13_model/programme/05_programme_report.py "${TAG:-}" && $PB/envs/pb_gpu_venv/bin/python -W ignore -u $PB/13_model/programme/04_annotate_programmes.py && echo "DONE all"
