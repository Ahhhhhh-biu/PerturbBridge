#!/bin/bash
#SBATCH --job-name=pb_bam_dl
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/bam_dl_%a.log
# One ENCODE genome BAM per task, md5-verified. List: 08_junction_psi/bam_list.tsv (url, dest, md5).
set -u
PB=/well/blagden/users/eic789/PerturbBridge; OUT=$PB/08_junction_psi/bam
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $PB/08_junction_psi/bam_list.tsv); url=$(echo "$line" | cut -f1); dest=$OUT/$(echo "$line" | cut -f2); md5=$(echo "$line" | cut -f3)
[ -s "$dest" ] && [ "$(md5sum "$dest" | cut -d' ' -f1)" = "$md5" ] && { echo "SKIP $dest"; exit 0; }
for try in 1 2 3 4 5 6; do
  curl -sS -L --retry 5 -m 7200 -o "$dest.part" "$url" && [ "$(md5sum "$dest.part" | cut -d' ' -f1)" = "$md5" ] && mv "$dest.part" "$dest" && echo "OK $dest" && exit 0
  echo "RETRY $try"; rm -f "$dest.part"; sleep $((try*30))
done
echo "FAIL $url"; exit 1
