#!/bin/bash
#SBATCH --job-name=pb_h4
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/h4_heya8.log
# Block 3 (H4 precancer cohorts) and block 4a (HeyA8 LARP1 KO validation), both need the dictionary (job 35037812)
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; export OMP_NUM_THREADS=8
$PB/envs/pb_gpu_venv/bin/python $PB/14_inverse/06_heya8_validation.py > $PB/logs/heya8_validation.txt 2>&1 && echo "DONE heya8 $(date)"
$PB/envs/pb_gpu_venv/bin/python $PB/14_inverse/05_h4_cohort_mapping.py > $PB/logs/h4_cohorts.txt 2>&1 && echo "DONE h4 $(date)"
