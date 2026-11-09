#!/bin/bash
#SBATCH --job-name=pb_ridge
#SBATCH --partition=gpu_a100_80gb,gpu_a100_40gb,gpu_rtx8000_48gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-19
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/ridge_%a.log
# ridge baseline: id 0-9 normal (cell*5+fold), 10-19 shuffled-RBP control
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; CL=(HepG2 K562); id=$SLURM_ARRAY_TASK_ID; j=$((id%10)); cl=${CL[$((j/5))]}; fold=$((j%5))
[ $id -lt 10 ] && M=normal || M=shuffle
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/benchmark/02_ridge_baseline.py $cl $fold --rbp_mode $M && echo "DONE ridge $cl $fold $M $(date)"
