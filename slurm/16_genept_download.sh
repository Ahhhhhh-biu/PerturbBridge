#!/bin/bash
#SBATCH --job-name=pb_genept
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/genept_download.log
# GenePT gene embeddings (Chen & Zou), Zenodo record 10833191, GenePT_emebdding_v2.zip (~574 MB; ada-002 and text-embedding-3-large pickles).
set -u
D=/well/blagden/users/eic789/PerturbBridge/06_features/rbp_annotation/genept; cd $D
curl -sS -L --fail --retry 5 -o GenePT_emebdding_v2.zip "https://zenodo.org/records/10833191/files/GenePT_emebdding_v2.zip?download=1" && unzip -o GenePT_emebdding_v2.zip && ls -la && echo "DONE genept" || echo "FAIL genept"
