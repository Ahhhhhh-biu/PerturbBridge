#!/bin/bash
#SBATCH --job-name=pb_gpuenv
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/gpu_env_build.log
# Builds envs/pb_gpu_venv: Python 3.11 venv with the CUDA build of torch (PyPI default wheel bundles the NVIDIA libraries)
# plus the analysis stack. Runs on a CPU node; the GPU check is a separate dependent job (21_gpu_env_test.sh).
set -u
PB=/well/blagden/users/eic789/PerturbBridge; export TMPDIR=$PB/envs/tmp; mkdir -p $TMPDIR
V=$PB/envs/pb_gpu_venv
[ -x $V/bin/python ] || /well/blagden/users/eic789/conda_envs/scrna_macro/bin/python -m venv $V
$V/bin/pip install --no-cache-dir --upgrade pip
$V/bin/pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cu128 && \
$V/bin/pip install --no-cache-dir numpy pandas pyarrow scipy scikit-learn statsmodels matplotlib pyyaml tqdm transformers safetensors pyfaidx pyBigWig && \
$V/bin/python -c "import torch, transformers; print('torch', torch.__version__, 'cuda build', torch.version.cuda)" && echo "DONE gpuenv" || echo "FAIL gpuenv"
