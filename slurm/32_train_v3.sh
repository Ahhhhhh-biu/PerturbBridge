#!/bin/bash
#SBATCH --job-name=pb_v3
#SBATCH --partition=gpu_a100_80gb,gpu_a100_40gb,gpu_rtx8000_48gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-49
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/v3_%a.log
# v3 = v2 recipe + --target clean (replicate-consistent, control-independent significant rows) + --no_psi (no run-specific control PSI input)
# id 0-29: cell*15 + fold*3 + tier (tiers A B C), tag _v3_s0 ; id 30-39: event-only control (tier A), tag _v3_ctrl_none_s0 ; id 40-49: shuffled-RBP control, tag _v3_ctrl_shuf_s0
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; export OMP_NUM_THREADS=4
CL=(HepG2 K562); TI=(A B C); id=$SLURM_ARRAY_TASK_ID
COMMON="--sel weighted --sigw 100 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20 --batch 2048 --seed 0 --target clean --no_psi"
if [ $id -lt 30 ]; then cl=${CL[$((id/15))]}; fold=$(( (id%15)/3 )); tier=${TI[$((id%3))]}; EXTRA="--tag _v3_s0"
elif [ $id -lt 40 ]; then j=$((id-30)); cl=${CL[$((j/5))]}; fold=$((j%5)); tier=A; EXTRA="--rbp_mode none --tag _v3_ctrl_none_s0"
else j=$((id-40)); cl=${CL[$((j/5))]}; fold=$((j%5)); tier=A; EXTRA="--rbp_mode shuffle --tag _v3_ctrl_shuf_s0"; fi
nvidia-smi --query-gpu=name --format=csv,noheader
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/01_train.py $cl $fold $tier $COMMON $EXTRA && echo "DONE v3 $cl $fold $tier $EXTRA $(date)"
