#!/bin/bash
#SBATCH --job-name=pb_cs
#SBATCH --partition=short
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/cohort_%x_%a.log
# One sample per task, streaming to keep disk use low: download its runs (40_sra_download_any.sh), fastp + STAR + rMATS prep (41_star_prep_any.sh),
# then delete the raw FASTQ of that sample once the prep is complete (public data, re-downloadable).
# Usage: sbatch --job-name=pb_cs_<GSE> --array=0-<n-1>%K --export=ALL,GSE=<GSE>,RL=<read length> 50_cohort_sample.sh
set -u
PB=/well/blagden/users/eic789/PerturbBridge; DIR=$PB/05_precancer_cohorts/$GSE; i=$SLURM_ARRAY_TASK_ID
line=$(awk -F'\t' -v i=$((i+2)) 'NR==i' $DIR/sample_table.tsv); hdr=$(head -1 $DIR/sample_table.tsv)
col() { echo "$hdr" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
gsm=$(echo "$line" | cut -f$(col gsm)); runs=$(echo "$line" | cut -f$(col runs) | tr ';' ' '); [ -n "$gsm" ] || exit 0
[ -f $DIR/rmats_prep/$gsm/.done ] && { echo "SKIP $gsm"; exit 0; }
for s in $runs; do idx=$(grep -o '^SRR[0-9]*' "$DIR/runinfo.csv" | sort -u | grep -n -x "$s" | cut -d: -f1); ( SLURM_ARRAY_TASK_ID=$((idx-1)); . $PB/slurm/40_sra_download_any.sh ); done
( SLURM_ARRAY_TASK_ID=$i; . $PB/slurm/41_star_prep_any.sh )
if [ -f $DIR/rmats_prep/$gsm/.done ]; then for s in $runs; do rm -f $DIR/fastq/${s}.fastq.gz $DIR/fastq/${s}_[12].fastq.gz $DIR/fastq/$s.ok; done; echo "CLEANED $gsm"; else echo "FAILED $gsm"; exit 1; fi
