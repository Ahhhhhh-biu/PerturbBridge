#!/bin/bash
#SBATCH --job-name=pb_bampost
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/bampost_%a.log
# rMATS-turbo post, one KD-vs-control pair per task (09_bam_layer/kd_control_pairs_all.tsv, header + 876 rows).
# Builds a per-pair tmp dir of symlinks to the four .rmats files from the prep step, lists the recorded BAM paths in --b1/--b2.
set -u
module load rMATS-turbo/4.3.0-gfbf-2023a
PB=/well/blagden/users/eic789/PerturbBridge; L=$PB/09_bam_layer; GTF=$PB/08_junction_psi/gencode.v29.gtf
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" $L/kd_control_pairs_all.tsv); [ -n "$line" ] || exit 0
kd=$(echo "$line" | cut -f1); ct=$(echo "$line" | cut -f2)
OUT=$L/rmats/${kd}_vs_${ct}; [ -f $OUT/.done ] && { echo "SKIP $kd"; exit 0; }
rm -rf $OUT; mkdir -p $OUT/tmp_post
b1=""; b2=""
for r in 1 2; do
  for e in $kd $ct; do
    P=$L/rmats_prep/${e}_rep$r; [ -f $P/.done ] || { echo "WAIT prep missing ${e}_rep$r"; exit 1; }
    for f in $P/tmp/*.rmats; do ln -s $f $OUT/tmp_post/${e}_rep${r}_$(basename $f); done
  done
  b1="$b1,$L/bam_tmp/${kd}_rep$r.bam"; b2="$b2,$L/bam_tmp/${ct}_rep$r.bam"
done
echo "${b1#,}" > $OUT/b1.txt; echo "${b2#,}" > $OUT/b2.txt
if rmats.py --b1 $OUT/b1.txt --b2 $OUT/b2.txt --gtf $GTF -t paired --readLength 100 --variable-read-length --allow-clipping --nthread 8 --task post --od $OUT --tmp $OUT/tmp_post > $OUT/post.log 2>&1 && [ -s $OUT/SE.MATS.JC.txt ]; then
  rm -rf $OUT/tmp_post $OUT/tmp; touch $OUT/.done; echo "DONE $kd vs $ct"
else echo "FAIL $kd vs $ct"; tail -5 $OUT/post.log; fi
