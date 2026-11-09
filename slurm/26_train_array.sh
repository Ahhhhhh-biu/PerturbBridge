#!/bin/bash
#SBATCH --job-name=pb_train
#SBATCH --partition=gpu_a100_80gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --time=03:50:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-29
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/train_%a.log
# task id -> (cell line, fold, tier): id = cell*15 + fold*3 + tier ; cells (HepG2, K562), folds 0-4, tiers (A, B, C)
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; export OMP_NUM_THREADS=4
CL=(HepG2 K562); TI=(A B C); id=$SLURM_ARRAY_TASK_ID; cl=${CL[$((id/15))]}; fold=$(( (id%15)/3 )); tier=${TI[$((id%3))]}
nvidia-smi --query-gpu=name --format=csv,noheader
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/01_train.py $cl $fold $tier --epochs 12 ${TRAIN_EXTRA:-} && echo "DONE train $cl $fold $tier $(date)"
