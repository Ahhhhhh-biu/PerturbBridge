#!/bin/bash
#SBATCH --job-name=pb_gputest
#SBATCH --partition=gpu_a100_80gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/gpu_test.log
# Checks whether the new gpu_blagden.prj account can actually be scheduled on a GPU node.
echo "node $(hostname) $(date)"; echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
module load CUDA 2>/dev/null && nvcc --version | tail -1
echo "DONE gputest"
