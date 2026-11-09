#!/bin/bash
#SBATCH --job-name=pb_seq
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=0-1
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/cnn_seq_%a.log
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; CL=(HepG2 K562)
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/benchmark/00_event_sequences.py ${CL[$SLURM_ARRAY_TASK_ID]}
