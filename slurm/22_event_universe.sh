#!/bin/bash
#SBATCH --job-name=pb_universe
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/event_universe_%j.log
# Full ENCODE event universe from all finished rMATS post runs (09_bam_layer/rmats/*/.done), then per-event sequence /
# conservation features and position-resolved eCLIP features. Output: 10_event_universe/encode_full/.
set -eu
module load BEDTools/2.31.0-GCC-12.3.0
PB=/well/blagden/users/eic789/PerturbBridge; U=$PB/10_event_universe/encode_full; PY=$PB/envs/pb_venv/bin/python
export TMPDIR=$PB/envs/tmp; mkdir -p $U
runs=$(for d in $PB/09_bam_layer/rmats/*/; do [ -f $d/.done ] && echo ${d%/}; done)
echo "runs: $(echo "$runs" | wc -l)  $(date)"
$PY $PB/10_event_universe/01_build_event_universe.py $U $runs; echo "universe built $(date)"
$PY $PB/10_event_universe/02_event_features.py $U; echo "features done $(date)"
$PY $PB/10_event_universe/03_eclip_event_features.py $U; echo "DONE universe $(date)"
