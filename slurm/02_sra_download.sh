#!/bin/bash
#SBATCH --job-name=pb_sra_dl
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/sra_dl_%x_%a.log
# NB: pass --job-name=pb_sra_<GSE> on the sbatch line so the two cohorts do not share log names.
# Usage: sbatch --array=0-N%6 --export=GSE=GSE274628 02_sra_download.sh
# Downloads one SRR per task: ENA fastq.gz (md5-verified) first, SRA Toolkit fallback. No home-directory cache.
set -u
PB=/well/blagden/users/eic789/PerturbBridge
DIR=$PB/05_precancer_cohorts/$GSE
OUT=$DIR/fastq; mkdir -p "$OUT"
SRR=$(grep -o '^SRR[0-9]*' "$DIR/runinfo.csv" | sort -u | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")
[ -n "$SRR" ] || { echo "no SRR for task $SLURM_ARRAY_TASK_ID"; exit 0; }
echo "$SRR start $(date)"
if ls "$OUT/${SRR}"_1.fastq.gz "$OUT/${SRR}"_2.fastq.gz >/dev/null 2>&1 && [ -f "$OUT/$SRR.ok" ]; then echo "SKIP $SRR"; exit 0; fi
rep=$(curl -sS -m 120 --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$SRR&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | tail -n +2)
ftp=$(echo "$rep" | cut -f2); md5s=$(echo "$rep" | cut -f3)
done_ena=0
if [ -n "$ftp" ] && [ "$(echo "$ftp" | tr ';' '\n' | wc -l)" -eq 2 ]; then
  i=0; allok=1
  for u in $(echo "$ftp" | tr ';' ' '); do
    i=$((i+1)); m=$(echo "$md5s" | cut -d';' -f$i); dest=$OUT/${SRR}_$i.fastq.gz
    if [ -s "$dest" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$m" ]; then echo "SKIP mate $i"; continue; fi
    for try in 1 2 3 4 5 6 7 8; do
      curl -sS -L --retry 5 -m 36000 -o "$dest.part" "https://$u" && [ "$(md5sum "$dest.part" | cut -d' ' -f1)" = "$m" ] && mv "$dest.part" "$dest" && break
      echo "RETRY $try $u"; rm -f "$dest.part"; sleep $((try*60))
    done
    [ -s "$dest" ] || allok=0
  done
  [ $allok -eq 1 ] && done_ena=1
fi
if [ $done_ena -eq 0 ]; then
  echo "ENA failed or incomplete for $SRR; falling back to SRA Toolkit"
  module load SRA-Toolkit/3.0.10-gompi-2023a
  export NCBI_SETTINGS=$PB/envs/ncbi_user-settings.mkfg
  TMP=$PB/05_precancer_cohorts/$GSE/tmp_$SRR; mkdir -p "$TMP"
  prefetch --max-size u -O "$TMP" "$SRR" && fasterq-dump --split-files --threads 4 --temp "$TMP" -O "$TMP" "$TMP/$SRR" \
    && gzip -c "$TMP/${SRR}_1.fastq" > "$OUT/${SRR}_1.fastq.gz" && gzip -c "$TMP/${SRR}_2.fastq" > "$OUT/${SRR}_2.fastq.gz" && rm -rf "$TMP"
fi
ls -s "$OUT/${SRR}"_1.fastq.gz "$OUT/${SRR}"_2.fastq.gz >/dev/null 2>&1 && touch "$OUT/$SRR.ok" && echo "OK $SRR $(date)" || echo "FAIL $SRR"
