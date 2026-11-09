#!/bin/bash
#SBATCH --job-name=pb_preptest
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/rmats_prep_post_test.log
# Test rMATS-turbo prep/post split on RBM25 shRNA (ENCSR610AEI) vs control (ENCSR603TCV):
# prep each BAM separately (via a symlink path, own tmp dir), delete the symlinks, link the .rmats files into one post tmp dir,
# run post with the BAM paths absent, and compare with the single-run ("both") result in 08_junction_psi/rmats.
set -u
module load rMATS-turbo/4.3.0-gfbf-2023a
PB=/well/blagden/users/eic789/PerturbBridge; D=$PB/08_junction_psi; T=$PB/09_bam_layer/prep_post_test; GTF=$D/gencode.v29.gtf
RM="-t paired --readLength 100 --variable-read-length --allow-clipping --nthread 8"
rm -rf $T; mkdir -p $T/links $T/tmp_post $T/post; cd $T
for b in ENCSR610AEI_rep1 ENCSR610AEI_rep2 ENCSR603TCV_rep1 ENCSR603TCV_rep2; do
  ln -s $D/bam/$b.bam $T/links/$b.bam; mkdir -p $T/prep_$b; echo "$T/links/$b.bam" > $T/prep_$b/b1.txt
  /usr/bin/time -f "prep $b %e s %M KB" rmats.py --b1 $T/prep_$b/b1.txt --gtf $GTF $RM --task prep --od $T/prep_$b --tmp $T/prep_$b/tmp || echo "FAIL prep $b"
  echo "prep $b files: $(ls $T/prep_$b/tmp) size $(du -sh $T/prep_$b/tmp | cut -f1); first line: $(head -c 200 $T/prep_$b/tmp/*_0.rmats | head -1)"
  for f in $T/prep_$b/tmp/*.rmats; do ln -s $f $T/tmp_post/${b}_$(basename $f); done
done
rm -f $T/links/*.bam; echo "symlinked BAM paths removed before post"
echo "$T/links/ENCSR610AEI_rep1.bam,$T/links/ENCSR610AEI_rep2.bam" > $T/b1.txt; echo "$T/links/ENCSR603TCV_rep1.bam,$T/links/ENCSR603TCV_rep2.bam" > $T/b2.txt
/usr/bin/time -f "post %e s %M KB" rmats.py --b1 $T/b1.txt --b2 $T/b2.txt --gtf $GTF $RM --task post --od $T/post --tmp $T/tmp_post || echo "FAIL post"
for t in SE A3SS A5SS MXE RI; do
  a=$T/post/$t.MATS.JC.txt; r=$D/rmats/ENCSR610AEI_vs_ENCSR603TCV/$t.MATS.JC.txt
  echo "$t post rows $(wc -l < $a) ref rows $(wc -l < $r); identical (coords+counts+PSI): $(cmp <(cut -f2-11,13-23 $a | sort) <(cut -f2-11,13-23 $r | sort) >/dev/null && echo yes || echo no)"
done
echo "TEST DONE"
