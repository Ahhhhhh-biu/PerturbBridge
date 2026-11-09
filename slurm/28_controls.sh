#!/bin/bash
#SBATCH --job-name=pb_ctrl
#SBATCH --partition=gpu_a100_80gb,gpu_a100_40gb,gpu_rtx8000_48gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-19
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/ctrl_%a.log
# RBP-specificity controls with the v2 training recipe, tier A, seed 0: id = mode*10 + cell*5 + fold
# mode 0 = none (event-only model, all RBP inputs zeroed), mode 1 = shuffle (each RBP gets another RBP's features)
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; export OMP_NUM_THREADS=4
MODE=(none shuffle); TAG=(_ctrl_none _ctrl_shuf); CL=(HepG2 K562); id=$SLURM_ARRAY_TASK_ID
m=$((id/10)); cl=${CL[$(((id%10)/5))]}; fold=$((id%5))
nvidia-smi --query-gpu=name --format=csv,noheader
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/01_train.py $cl $fold A --sel weighted --sigw 100 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20 --batch 2048 --seed 0 --rbp_mode ${MODE[$m]} --tag ${TAG[$m]}_s0 && echo "DONE ctrl ${MODE[$m]} $cl $fold $(date)"
