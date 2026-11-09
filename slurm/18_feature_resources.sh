#!/bin/bash
#SBATCH --job-name=pb_featres
#SBATCH --partition=short
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/feature_resources.log
# Downloads event-encoder resources: MaxEntScan (Yeo & Burge 2004) scripts, UCSC hg38 phastCons100way bigWig, and bigWigAverageOverBed.
set -u
R=/well/blagden/users/eic789/PerturbBridge/06_features/resources; cd $R
get(){ for i in 1 2 3 4 5 6; do curl -sSL --retry 3 -o "$2" "$1" && return 0; sleep 30; done; return 1; }
get http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz maxentscan.tar.gz && mkdir -p maxentscan && tar -xzf maxentscan.tar.gz -C maxentscan && echo "DONE maxentscan $(ls maxentscan/fordownload | head -3 | tr '\n' ' ')" || echo "FAIL maxentscan"
get https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw hg38.phastCons100way.bw && echo "DONE phastCons $(du -h hg38.phastCons100way.bw | cut -f1)" || echo "FAIL phastCons"
get https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed bigWigAverageOverBed && chmod +x bigWigAverageOverBed && echo "DONE bigWigAverageOverBed" || echo "FAIL bigWigAverageOverBed"
