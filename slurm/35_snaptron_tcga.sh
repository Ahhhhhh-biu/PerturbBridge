#!/bin/bash
#SBATCH --job-name=pb_snap
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --array=0-47
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/snaptron_%a.log
# TCGA per-sample PSI for the model event universe via Snaptron (14_inverse/01_snaptron_tcga_psi.py), 48 shards
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp
$PB/envs/pb_gpu_venv/bin/python $PB/14_inverse/01_snaptron_tcga_psi.py $SLURM_ARRAY_TASK_ID ${NTASK:-$SLURM_ARRAY_TASK_COUNT}
