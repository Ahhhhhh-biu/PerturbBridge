#!/bin/bash
#SBATCH --job-name=pb_cnn
#SBATCH --partition=gpu_a100_80gb,gpu_a100_40gb,gpu_rtx8000_48gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --time=02:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-19
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/cnn_%a.log
# v4 CNN benchmark: id 0-9 = cell*5 + fold, tier A, tag _v4cnn_s0 ; id 10-19 = shuffled-RBP control, tag _v4cnn_ctrl_shuf_s0. Same recipe as v3.
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; export OMP_NUM_THREADS=4
CL=(HepG2 K562); id=$SLURM_ARRAY_TASK_ID; j=$((id%10)); cl=${CL[$((j/5))]}; fold=$((j%5))
COMMON="--sel weighted --sigw 100 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20 --batch 2048 --seed 0 --target clean --no_psi"
if [ $id -lt 10 ]; then EXTRA="--tag _v4cnn_s0"; else EXTRA="--rbp_mode shuffle --tag _v4cnn_ctrl_shuf_s0"; fi
nvidia-smi --query-gpu=name --format=csv,noheader
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/benchmark/01_train_cnn.py $cl $fold A $COMMON $EXTRA && echo "DONE cnn $cl $fold $EXTRA $(date)"
