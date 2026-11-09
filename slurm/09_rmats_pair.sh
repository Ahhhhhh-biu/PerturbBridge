#!/bin/bash
#SBATCH --job-name=pb_rmats
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/rmats_%a.log
# rMATS-turbo, one KD-vs-control pair per task (08_junction_psi/kd_control_pairs.tsv), GENCODE v29, paired-end, variable read length.
set -u
module load rMATS-turbo/4.3.0-gfbf-2023a
PB=/well/blagden/users/eic789/PerturbBridge; D=$PB/08_junction_psi
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" $D/kd_control_pairs.tsv); kd=$(echo "$line" | cut -f1); ct=$(echo "$line" | cut -f2)
[ -n "$kd" ] || exit 0
OUT=$D/rmats/${kd}_vs_${ct}; mkdir -p $OUT/tmp
ls $D/bam/${kd}_rep*.bam | paste -sd, > $OUT/b1.txt; ls $D/bam/${ct}_rep*.bam | paste -sd, > $OUT/b2.txt
[ -s $D/gencode.v29.gtf ] || zcat $PB/03_suppa_psi/gencode.v29.primary_assembly.annotation.gtf.gz > $D/gencode.v29.gtf
echo "$kd vs $ct: b1=$(cat $OUT/b1.txt) b2=$(cat $OUT/b2.txt)"
rmats.py --b1 $OUT/b1.txt --b2 $OUT/b2.txt --gtf $D/gencode.v29.gtf -t paired --readLength 100 --variable-read-length --allow-clipping --nthread 8 --od $OUT --tmp $OUT/tmp && rm -rf $OUT/tmp && echo "DONE $kd" || echo "FAIL $kd"
