#!/bin/bash
#SBATCH --job-name=pb_eclip_feat
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/eclip_features.log
module load BEDTools/2.31.0-GCC-12.3.0
PB=/well/blagden/users/eic789/PerturbBridge
$PB/envs/pb_venv/bin/python $PB/06_features/01_eclip_event_features.py
