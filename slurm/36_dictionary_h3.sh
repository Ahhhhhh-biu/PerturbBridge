#!/bin/bash
#SBATCH --job-name=pb_dict
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/dictionary_h3.log
# Block 1 (dictionary) then block 2a (H3 retrieval inside ENCODE)
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; export OMP_NUM_THREADS=8
$PB/envs/pb_gpu_venv/bin/python $PB/14_inverse/02_build_dictionary.py && echo "DONE dictionary $(date)" && \
$PB/envs/pb_gpu_venv/bin/python $PB/14_inverse/03_h3_encode_retrieval.py > $PB/logs/h3_encode.txt 2>&1 && echo "DONE h3_encode $(date)"
