#!/bin/bash
#SBATCH --job-name=pb_atlas_esm2
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/atlas_esm2.txt
# ESM2 embeddings for atlas proteins outside the ENCODE set.
PB=/well/blagden/users/eic789/PerturbBridge; export HF_HOME=$PB/envs/hf_cache
$PB/envs/pb_gpu_venv/bin/python -W ignore -u $PB/13_model/programme/08_atlas_esm2.py && echo "DONE atlas esm2"
