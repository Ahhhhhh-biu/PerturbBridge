#!/bin/bash
#SBATCH --job-name=pb_sra_dl
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/sra_dl_%x_%a.log
# Single- or paired-end variant of 02_sra_download.sh. Usage: sbatch --job-name=pb_sra_<GSE> --array=0-N%10 --export=ALL,GSE=<GSE> 40_sra_download_any.sh
# One SRR per task: ENA fastq.gz (md5-verified; 1 file for single-end, 2 for paired) with SRA Toolkit fallback. Writes <SRR>.ok on success.
set -u
PB=/well/blagden/users/eic789/PerturbBridge; DIR=$PB/05_precancer_cohorts/$GSE; OUT=$DIR/fastq; mkdir -p "$OUT"
SRR=$(grep -o '^SRR[0-9]*' "$DIR/runinfo.csv" | sort -u | sed -n "$((SLURM_ARRAY_TASK_ID+1))p"); [ -n "$SRR" ] || { echo "no SRR for task $SLURM_ARRAY_TASK_ID"; exit 0; }
[ -f "$OUT/$SRR.ok" ] && { echo "SKIP $SRR"; exit 0; }; echo "$SRR start $(date)"
rep=$(curl -sS -m 120 --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRR&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | tail -n +2)
ftp=$(echo "$rep" | cut -f2); md5s=$(echo "$rep" | cut -f3); nf=$(echo "$ftp" | tr ';' '\n' | grep -c .); ok=0
if [ "$nf" -ge 1 ]; then
  i=0; allok=1
  for u in $(echo "$ftp" | tr ';' ' '); do
    i=$((i+1)); m=$(echo "$md5s" | cut -d';' -f$i); if [ "$nf" -eq 1 ]; then dest=$OUT/${SRR}.fastq.gz; else dest=$OUT/${SRR}_$i.fastq.gz; fi
    if [ -s "$dest" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$m" ]; then echo "SKIP file $i"; continue; fi
    for try in 1 2 3 4 5 6; do
      curl -sS -L --retry 5 -m 36000 -o "$dest.part" "https://$u" && [ "$(md5sum "$dest.part" | cut -d' ' -f1)" = "$m" ] && mv "$dest.part" "$dest" && break
      echo "RETRY $try $u"; rm -f "$dest.part"; sleep $((try*60))
    done
    [ -s "$dest" ] || allok=0
  done
  [ $allok -eq 1 ] && ok=1
fi
if [ $ok -eq 0 ]; then
  echo "ENA failed for $SRR; SRA Toolkit fallback"; module load SRA-Toolkit/3.0.10-gompi-2023a; export NCBI_SETTINGS=$PB/envs/ncbi_user-settings.mkfg
  TMP=$DIR/tmp_$SRR; mkdir -p "$TMP"
  if prefetch --max-size u -O "$TMP" "$SRR" && fasterq-dump --split-files --threads 4 --temp "$TMP" -O "$TMP" "$TMP/$SRR"; then
    if [ -f "$TMP/${SRR}_2.fastq" ]; then gzip -c "$TMP/${SRR}_1.fastq" > "$OUT/${SRR}_1.fastq.gz" && gzip -c "$TMP/${SRR}_2.fastq" > "$OUT/${SRR}_2.fastq.gz" && ok=1
    else f=$(ls "$TMP"/${SRR}*.fastq | head -1); gzip -c "$f" > "$OUT/${SRR}.fastq.gz" && ok=1; fi
  fi; rm -rf "$TMP"
fi
[ $ok -eq 1 ] && touch "$OUT/$SRR.ok" && echo "OK $SRR $(date)" || echo "FAIL $SRR"
