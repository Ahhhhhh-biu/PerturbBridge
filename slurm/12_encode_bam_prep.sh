#!/bin/bash
#SBATCH --job-name=pb_bamprep
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/bamprep_%a.log
# Streaming ENCODE BAM layer, one BAM per task (09_bam_layer/bam_manifest_needed.tsv, header + 1956 rows):
# download (or reuse a BAM already in 08_junction_psi/bam) -> md5 -> rMATS-turbo prep (GENCODE v29) -> delete the downloaded BAM.
# Keeps only 09_bam_layer/rmats_prep/<dest>/tmp/*.rmats and a .done marker. The BAM path recorded inside .rmats is
# 09_bam_layer/bam_tmp/<dest>.bam, which is what the post step must list in --b1/--b2.
set -u
module load rMATS-turbo/4.3.0-gfbf-2023a
PB=/well/blagden/users/eic789/PerturbBridge; L=$PB/09_bam_layer; M=$L/bam_manifest_needed.tsv; GTF=$PB/08_junction_psi/gencode.v29.gtf
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" $M); [ -n "$line" ] || exit 0
url=$(echo "$line" | cut -f1); dest=$(echo "$line" | cut -f2); md5=$(echo "$line" | cut -f3); name=${dest%.bam}
P=$L/rmats_prep/$name; B=$L/bam_tmp/$dest
[ -f $P/.done ] && { echo "SKIP $name"; exit 0; }
mkdir -p $P $L/bam_tmp
downloaded=0
if [ -s $PB/08_junction_psi/bam/$dest ]; then ln -sf $PB/08_junction_psi/bam/$dest $B; echo "reuse $dest"
else
  ok=0
  for try in 1 2 3 4 5 6; do
    curl -sS -L --fail --retry 3 -C - -o $B.part "$url" && [ "$(md5sum $B.part | cut -d' ' -f1)" = "$md5" ] && { mv $B.part $B; ok=1; break; }
    echo "retry $try $name"; sleep $((try*30))
  done
  [ $ok -eq 1 ] || { echo "FAIL download $name"; rm -f $B.part; exit 1; }
  downloaded=1
fi
echo "$B" > $P/b1.txt; rm -rf $P/tmp
if rmats.py --b1 $P/b1.txt --gtf $GTF -t paired --readLength 100 --variable-read-length --allow-clipping --nthread 4 --task prep --od $P --tmp $P/tmp > $P/prep.log 2>&1 \
   && ls $P/tmp/*.rmats >/dev/null 2>&1; then
  rm -rf $P/tmp/*bam* $P/tmp/*.txt 2>/dev/null; touch $P/.done; echo "DONE $name $(date) rmats $(du -sh $P/tmp | cut -f1)"
else
  echo "FAIL prep $name"; tail -5 $P/prep.log
fi
rm -f $B
