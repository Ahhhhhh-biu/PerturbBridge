#!/bin/bash
# Usage: bash 02b_resubmit_failed.sh GSE274628 [comma-separated task indices]; without indices, all runs lacking .ok are resubmitted.
GSE=$1; PB=/well/blagden/users/eic789/PerturbBridge; DIR=$PB/05_precancer_cohorts/$GSE
if [ -n "$2" ]; then idx=$2; else
  idx=$(grep -o '^SRR[0-9]*' $DIR/runinfo.csv | sort -u | awk -v d="$DIR/fastq" '{ if (system("test -f " d "/" $1 ".ok")!=0) print NR-1 }' | paste -sd,)
fi
[ -z "$idx" ] && { echo "nothing to resubmit for $GSE"; exit 0; }
for i in $(echo $idx | tr ',' ' '); do rm -rf $DIR/tmp_$(grep -o '^SRR[0-9]*' $DIR/runinfo.csv | sort -u | sed -n "$((i+1))p"); done
j=$(sbatch --parsable --job-name=pb_sra_$GSE --array=$idx%6 --export=ALL,GSE=$GSE $PB/slurm/02_sra_download.sh)
echo "$j ${GSE}_fastq resubmit indices $idx" | tee -a $PB/logs/JOBS.txt
