#!/bin/bash
#SBATCH --job-name=pb_starprep
#SBATCH --partition=short
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/starprep_%x_%a.log
# Usage: sbatch --job-name=pb_starprep_<GSE> --array=0-<n-1>%8 --export=ALL,GSE=<GSE>,RL=<read length> 14_star_prep_cohort.sh
# Per sample: fastp adapter/polyG trimming (paired, auto-detect, min length 36) -> STAR (GRCh38 + GENCODE v47 index) -> rMATS-turbo prep (v47 GTF) -> delete trimmed FASTQ and BAM. Keeps SJ.out.tab, ReadsPerGene,
# Log.final.out and rmats_prep/<gsm>/tmp/*.rmats. Disk-safe variant of 07_star_align.sh for the 94-sample fallopian cohort.
set -u
module load fastp/0.23.4-GCC-12.3.0 STAR/2.7.11b-GCC-12.3.0 SAMtools/1.18-GCC-12.3.0 rMATS-turbo/4.3.0-gfbf-2023a
PB=/well/blagden/users/eic789/PerturbBridge; DIR=$PB/05_precancer_cohorts/$GSE; FQ=$DIR/fastq; OUT=$DIR/star; mkdir -p $OUT
IDX=/well/blagden/users/eic789/HeyA8_lnc_antigen/02_reference/STAR_GRCh38_v47; GTF=/well/blagden/users/eic789/LncRNASPLICE/gencode.v47.primary_assembly.annotation.gtf
line=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID+2)) 'NR==i' $DIR/sample_table.tsv); hdr=$(head -1 $DIR/sample_table.tsv)
col() { echo "$hdr" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
gsm=$(echo "$line" | cut -f$(col gsm)); runs=$(echo "$line" | cut -f$(col runs) | tr ';' ' ')
[ -n "$gsm" ] || exit 0
P=$DIR/rmats_prep/$gsm; [ -f $P/.done ] && { echo "SKIP $gsm"; exit 0; }
r1=""; r2=""
for s in $runs; do [ -f $FQ/$s.ok ] || { echo "INCOMPLETE $gsm: run $s"; exit 1; }; r1="$r1,$FQ/${s}_1.fastq.gz"; r2="$r2,$FQ/${s}_2.fastq.gz"; done
BAM=$OUT/$gsm.Aligned.sortedByCoord.out.bam; TR=$DIR/trimmed_tmp/$gsm; mkdir -p $TR $DIR/fastp
if [ ! -f $OUT/$gsm.done ]; then
  echo "$gsm fastp start $(date)"
  fastp -w 8 -i <(zcat ${r1//,/ }) -I <(zcat ${r2//,/ }) -o $TR/R1.fq.gz -O $TR/R2.fq.gz --detect_adapter_for_pe --trim_poly_g --length_required 36 \
    -j $DIR/fastp/$gsm.fastp.json -h $DIR/fastp/$gsm.fastp.html 2> $DIR/fastp/$gsm.fastp.log || { echo "FAIL fastp $gsm"; exit 1; }
  echo "$gsm STAR start $(date)"
  STAR --runThreadN 8 --genomeDir $IDX --readFilesIn $TR/R1.fq.gz $TR/R2.fq.gz --readFilesCommand zcat \
    --outFileNamePrefix $OUT/$gsm. --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
    --outSAMattributes NH HI AS nM NM XS --twopassMode Basic --quantMode GeneCounts \
    --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --outSAMattrRGline ID:$gsm SM:$gsm \
    --limitBAMsortRAM 30000000000 --outTmpDir $OUT/tmp_$gsm || { echo "FAIL STAR $gsm"; exit 1; }
  rm -rf $OUT/tmp_$gsm $OUT/$gsm._STARgenome $OUT/$gsm._STARpass1 $TR; samtools index $BAM && touch $OUT/$gsm.done
fi
mkdir -p $P; echo "$BAM" > $P/b1.txt; rm -rf $P/tmp
if rmats.py --b1 $P/b1.txt --gtf $GTF -t paired --readLength $RL --variable-read-length --allow-clipping --nthread 8 --task prep --od $P --tmp $P/tmp > $P/prep.log 2>&1 && ls $P/tmp/*.rmats >/dev/null 2>&1; then
  touch $P/.done; rm -f $BAM $BAM.bai; echo "DONE $gsm $(date) uniq $(grep -m1 'Uniquely mapped reads %' $OUT/$gsm.Log.final.out | cut -f2)"
else echo "FAIL prep $gsm"; tail -3 $P/prep.log; fi
