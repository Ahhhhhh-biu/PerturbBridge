#!/bin/bash
#SBATCH --job-name=pb_baselines
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=0-3
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/baselines_%a.log
# tasks: 0 HepG2 raw dPSI, 1 K562 raw, 2 HepG2 control-centred, 3 K562 control-centred
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; export OMP_NUM_THREADS=8
CL=(HepG2 K562 HepG2 K562); LAB=(dpsi dpsi dpsi_cc dpsi_cc)
$PB/envs/pb_venv/bin/python $PB/12_baselines/02_baselines.py ${CL[$SLURM_ARRAY_TASK_ID]} ${LAB[$SLURM_ARRAY_TASK_ID]} && echo "DONE baselines ${CL[$SLURM_ARRAY_TASK_ID]} ${LAB[$SLURM_ARRAY_TASK_ID]} $(date)"
