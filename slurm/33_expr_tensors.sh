#!/bin/bash
#SBATCH --job-name=pb_exprT
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-1
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/expr_tensors_%a.log
# gene-expression layer tensors + measurement references per cell line (13_model/10_build_expr_tensors.py)
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; module load BEDTools/2.31.0-GCC-12.3.0
CL=(HepG2 K562); cl=${CL[$SLURM_ARRAY_TASK_ID]}
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/10_build_expr_tensors.py $cl && echo "DONE expr tensors $cl $(date)"
