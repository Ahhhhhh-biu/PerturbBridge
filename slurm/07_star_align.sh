#!/bin/bash
#SBATCH --job-name=pb_star
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/star_%x_%a.log
# Usage: sbatch --job-name=pb_star_<GSE> --array=0-<n-1>%8 --export=ALL,GSE=<GSE> 07_star_align.sh
# One sample per task (all its SRR runs given to STAR as comma-separated lists). Skips samples whose runs are not all downloaded.
set -u
module load STAR/2.7.11b-GCC-12.3.0 SAMtools/1.18-GCC-12.3.0
PB=/well/blagden/users/eic789/PerturbBridge
DIR=$PB/05_precancer_cohorts/$GSE; FQ=$DIR/fastq; OUT=$DIR/star; mkdir -p $OUT
IDX=/well/blagden/users/eic789/HeyA8_lnc_antigen/02_reference/STAR_GRCh38_v47
line=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID+2)) 'NR==i' $DIR/sample_table.tsv)
hdr=$(head -1 $DIR/sample_table.tsv)
col() { echo "$hdr" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
gsm=$(echo "$line" | cut -f$(col gsm)); runs=$(echo "$line" | cut -f$(col runs) | tr ';' ' ')
[ -n "$gsm" ] || { echo "no sample for task $SLURM_ARRAY_TASK_ID"; exit 0; }
[ -s $OUT/$gsm.Aligned.sortedByCoord.out.bam ] && [ -f $OUT/$gsm.done ] && { echo "SKIP $gsm"; exit 0; }
r1=""; r2=""
for s in $runs; do [ -f $FQ/$s.ok ] || { echo "INCOMPLETE $gsm: run $s not downloaded"; exit 0; }; r1="$r1,$FQ/${s}_1.fastq.gz"; r2="$r2,$FQ/${s}_2.fastq.gz"; done
r1=${r1#,}; r2=${r2#,}
echo "$gsm runs: $runs  start $(date)"
STAR --runThreadN 8 --genomeDir $IDX --readFilesIn $r1 $r2 --readFilesCommand zcat \
  --outFileNamePrefix $OUT/$gsm. --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
  --outSAMattributes NH HI AS nM NM XS --twopassMode Basic --quantMode GeneCounts \
  --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --outSAMattrRGline ID:$gsm SM:$gsm \
  --limitBAMsortRAM 30000000000 --outTmpDir $OUT/tmp_$gsm
rc=$?
if [ $rc -eq 0 ]; then samtools index $OUT/$gsm.Aligned.sortedByCoord.out.bam && touch $OUT/$gsm.done && echo "DONE $gsm $(date)"; else echo "FAIL $gsm rc=$rc"; fi
rm -rf $OUT/tmp_$gsm $OUT/$gsm._STARgenome $OUT/$gsm._STARpass1
