#!/bin/bash
#SBATCH --job-name=pb_ftpost
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/rmats_post_GSE274628.log
# rMATS-turbo post for the fallopian cohort GSE274628: b1 = BRCA1/2 carriers (41), b2 = non-carriers (53).
# Uses the per-sample prep outputs (14_star_prep_cohort.sh), symlinked into tmp_post; BAMs were deleted after prep.
set -u
module load rMATS-turbo/4.3.0-gfbf-2023a
PB=/well/blagden/users/eic789/PerturbBridge; OUT=$PB/05_precancer_cohorts/GSE274628/rmats_carrier_vs_noncarrier
GTF=/well/blagden/users/eic789/LncRNASPLICE/gencode.v47.primary_assembly.annotation.gtf
echo "start $(date)"
if rmats.py --b1 $OUT/b1.txt --b2 $OUT/b2.txt --gtf $GTF -t paired --readLength 101 --variable-read-length --allow-clipping --nthread 16 --task post --od $OUT --tmp $OUT/tmp_post > $OUT/post.log 2>&1 && [ -s $OUT/SE.MATS.JC.txt ]; then
  rm -rf $OUT/tmp_post; touch $OUT/.done; echo "DONE $(date)"
else echo "FAIL"; tail -5 $OUT/post.log; fi
