#!/bin/bash
#SBATCH --job-name=pb_esm2
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/esm2_embeddings.log
# Installs torch (CPU wheels) + transformers into pb_venv (TMPDIR under /well), downloads ESM2 650M via the cluster proxy, embeds all RBP sequences.
set -u
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; mkdir -p $TMPDIR $PB/envs/hf_cache; export HF_HOME=$PB/envs/hf_cache
$PB/envs/pb_venv/bin/python -c "import torch, transformers" 2>/dev/null || $PB/envs/pb_venv/bin/pip install --no-cache-dir --index-url https://download.pytorch.org/whl/cpu torch && $PB/envs/pb_venv/bin/pip install --no-cache-dir transformers safetensors
$PB/envs/pb_venv/bin/python -c "import torch, transformers; print('torch', torch.__version__, 'transformers', transformers.__version__)" || { echo "FAIL install"; exit 1; }
$PB/envs/pb_venv/bin/python $PB/06_features/rbp_annotation/02_esm2_embeddings.py && echo "DONE esm2" || echo "FAIL esm2"
