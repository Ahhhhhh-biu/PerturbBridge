#!/bin/bash
#SBATCH --job-name=pb_suppa_events
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/suppa_events.log
# Generate SUPPA2 local AS events (SE, SS, MX, RI, FL) and transcript-level ioi from GENCODE v29 (ENCODE's annotation).
set -eu
PB=/well/blagden/users/eic789/PerturbBridge
PY=$PB/envs/pb_venv/bin/python
D=$PB/03_suppa_psi
zcat $D/gencode.v29.primary_assembly.annotation.gtf.gz > $D/gencode.v29.gtf
$PY $PB/envs/SUPPA/suppa.py generateEvents -i $D/gencode.v29.gtf -o $D/events/gencode_v29 -f ioe -e SE SS MX RI FL
$PY $PB/envs/SUPPA/suppa.py generateEvents -i $D/gencode.v29.gtf -o $D/events/gencode_v29 -f ioi
cat $D/events/gencode_v29_*_strict.ioe | awk 'NR==1 || !/^seqname/' > $D/events/gencode_v29_all_strict.ioe
wc -l $D/events/*.ioe $D/events/*.ioi
rm $D/gencode.v29.gtf
