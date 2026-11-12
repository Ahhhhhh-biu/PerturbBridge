#!/bin/bash
#SBATCH --job-name=pb_finish
#SBATCH --partition=long
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/cohort_finish_%x.log
# Controller for one cohort: wait for the running sample array, resubmit samples without a finished prep (up to 3 rounds),
# then run rMATS post over all samples and the inverse-mapping analysis. Stops if samples are still missing after the retries.
# Usage: sbatch --job-name=pb_finish_<GSE> --export=ALL,GSE=<GSE>,RL=<len>,LAYOUT=paired|single,K=<concurrency>,WAITFOR=<array job id> 54_cohort_finish.sh
PB=/well/blagden/users/eic789/PerturbBridge; DIR=$PB/05_precancer_cohorts/$GSE; n=$(($(wc -l < $DIR/sample_table.tsv)-1)); cd $PB
while squeue -h -j $WAITFOR 2>/dev/null | grep -q .; do sleep 300; done
missing() { m=""; for i in $(seq 0 $((n-1))); do g=$(awk -F'\t' -v i=$((i+2)) 'NR==i{print $1}' $DIR/sample_table.tsv); [ -f $DIR/rmats_prep/$g/.done ] || m="$m,$i"; done; echo "${m#,}"; }
for round in 1 2 3; do m=$(missing); [ -z "$m" ] && break; echo "round $round: resubmitting $(echo $m | tr ',' '\n' | wc -l) samples $(date)"; sbatch --wait --job-name=pb_cs_$GSE --array=$m%$K --export=ALL,GSE=$GSE,RL=$RL $PB/slurm/50_cohort_sample.sh; done
m=$(missing); [ -z "$m" ] || { echo "STOP: samples still missing after retries: $m"; exit 1; }
echo "all $n samples prepared $(date)"; rm -rf $DIR/rmats_all/.done
sbatch --wait --time=24:00:00 --job-name=pb_rmpost_$GSE --export=ALL,GSE=$GSE,NAME=rmats_all,RL=$RL,LAYOUT=$LAYOUT $PB/slurm/45_rmats_post_any.sh
sbatch --wait --job-name=pb_h4b2_$GSE --export=ALL,GSE=$GSE $PB/slurm/53_h4_batch2.sh; echo "FINISHED $GSE $(date)"; tail -3 $PB/logs/h4_batch2_$GSE.txt
