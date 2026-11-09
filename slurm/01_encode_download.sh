#!/bin/bash
#SBATCH --job-name=pb_encode_dl
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --array=0-67%10
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/encode_dl_%a.log
# Download ENCODE quantification TSVs + eCLIP peaks listed in download_manifest.tsv (100 files per array task), verify md5.
set -u
PB=/well/blagden/users/eic789/PerturbBridge
MAN=$PB/01_encode_metadata/download_list.tsv   # url<TAB>dest<TAB>md5, no header
OUT=$PB/02_encode_download
CHUNK=100
start=$(( SLURM_ARRAY_TASK_ID * CHUNK + 1 ))
end=$(( start + CHUNK - 1 ))
sed -n "${start},${end}p" "$MAN" | while IFS=$'\t' read -r url d md5; do
  dest=$OUT/$d
  if [ -s "$dest" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$md5" ]; then echo "SKIP $dest"; continue; fi
  for try in 1 2 3 4 5 6; do
    curl -sS -L --retry 5 -m 1800 -o "$dest.part" "$url" && [ "$(md5sum "$dest.part" | cut -d' ' -f1)" = "$md5" ] && mv "$dest.part" "$dest" && echo "OK $dest" && break
    echo "RETRY $try $url"; sleep $((try*20))
  done
  [ -s "$dest" ] || echo "FAIL $url"
done
echo "task $SLURM_ARRAY_TASK_ID done $(date)"
