#!/bin/bash
#SBATCH --job-name=pb_expr
#SBATCH --partition=gpu_a100_80gb,gpu_a100_40gb,gpu_rtx8000_48gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-39
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/expr_%a.log
# gene-expression output layer: id = cell*20 + variant*5 + fold ; variants 0 = tier A (_expr_s0), 1 = tier B (_expr_s0), 2 = gene-only control (_expr_ctrl_none_s0), 3 = shuffled-RBP control (_expr_ctrl_shuf_s0)
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; export OMP_NUM_THREADS=4
CL=(HepG2 K562); id=$SLURM_ARRAY_TASK_ID; cl=${CL[$((id/20))]}; v=$(( (id%20)/5 )); fold=$((id%5))
COMMON="--epochs 20 --batch 4096 --lr 3e-4 --sigw 10 --rpca 32 --rdrop 0.3 --cos --seed 0"
case $v in 0) tier=A; EXTRA="--tag _expr_s0";; 1) tier=B; EXTRA="--tag _expr_s0";; 2) tier=A; EXTRA="--rbp_mode none --tag _expr_ctrl_none_s0";; 3) tier=A; EXTRA="--rbp_mode shuffle --tag _expr_ctrl_shuf_s0";; esac
nvidia-smi --query-gpu=name --format=csv,noheader
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/11_train_expr.py $cl $fold $tier $COMMON $EXTRA && echo "DONE expr $cl $fold $tier $EXTRA $(date)"
