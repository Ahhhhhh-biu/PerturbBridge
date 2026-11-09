#!/bin/bash
#SBATCH --job-name=pb_rmpost
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/rmats_post_%x.log
# rMATS-turbo post for any cohort prepared with 41_star_prep_any.sh. Usage:
#   sbatch --job-name=pb_rmpost_<GSE> --export=ALL,GSE=<GSE>,NAME=<out dir name under the cohort>,RL=<read length>,LAYOUT=single|paired 45_rmats_post_any.sh
# <cohort>/<NAME>/b1_samples.txt and b2_samples.txt (one GSM per line) must exist. The per-sample .rmats files from rmats_prep/<gsm>/tmp
# are symlinked into tmp_post; post matches them by the BAM path recorded in each .rmats file, so the (deleted) BAMs are not needed.
set -u
module load rMATS-turbo/4.3.0-gfbf-2023a
PB=/well/blagden/users/eic789/PerturbBridge; DIR=$PB/05_precancer_cohorts/$GSE; OUT=$DIR/$NAME
GTF=/well/blagden/users/eic789/LncRNASPLICE/gencode.v47.primary_assembly.annotation.gtf
rm -rf $OUT/tmp_post; mkdir -p $OUT/tmp_post; missing=0
for grp in b1 b2; do
  list=""
  for g in $(cat $OUT/${grp}_samples.txt); do
    [ -f $DIR/rmats_prep/$g/.done ] || { echo "MISSING prep $g"; missing=$((missing+1)); continue; }
    for f in $DIR/rmats_prep/$g/tmp/*.rmats; do ln -s $f $OUT/tmp_post/${g}_$(basename $f); done
    list="$list,$DIR/star/$g.Aligned.sortedByCoord.out.bam"
  done
  echo "${list#,}" > $OUT/$grp.txt
done
echo "start $(date) b1 $(tr ',' '\n' < $OUT/b1.txt | wc -l) b2 $(tr ',' '\n' < $OUT/b2.txt | wc -l) missing $missing"
if rmats.py --b1 $OUT/b1.txt --b2 $OUT/b2.txt --gtf $GTF -t $LAYOUT --readLength $RL --variable-read-length --allow-clipping --nthread 16 --task post --od $OUT --tmp $OUT/tmp_post > $OUT/post.log 2>&1 && [ -s $OUT/SE.MATS.JC.txt ]; then
  rm -rf $OUT/tmp_post; touch $OUT/.done; echo "DONE $GSE $NAME $(date) (missing prep: $missing)"
else echo "FAIL"; tail -5 $OUT/post.log; fi
