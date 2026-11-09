#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
# Usage: sbatch --job-name=pb_rmats_<GSE> --output=.../logs/rmats_cohort_<GSE>.log --export=ALL,OUT=<dir with b1.txt b2.txt>,RL=<read length> 11_rmats_cohort.sh
# Group-level rMATS-turbo (both tasks) on STAR GRCh38 + GENCODE v47 alignments: b1 = case group, b2 = control group.
set -u
module load rMATS-turbo/4.3.0-gfbf-2023a
GTF=/well/blagden/users/eic789/LncRNASPLICE/gencode.v47.primary_assembly.annotation.gtf
echo "start $(date) OUT=$OUT RL=$RL"
rmats.py --b1 $OUT/b1.txt --b2 $OUT/b2.txt --gtf $GTF -t paired --readLength $RL --variable-read-length --allow-clipping --nthread 16 --od $OUT --tmp $OUT/tmp && rm -rf $OUT/tmp && echo "DONE $(date)" || echo "FAIL"
