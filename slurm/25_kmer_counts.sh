#!/bin/bash
#SBATCH --job-name=pb_kmer
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/kmer_counts_%j.log
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp
$PB/envs/pb_venv/bin/python $PB/10_event_universe/06_kmer_counts.py $PB/10_event_universe/encode_full && echo "DONE kmer $(date)"
