#!/bin/bash
#SBATCH --job-name=pb_gpuenvtest
#SBATCH --partition=gpu_a100_80gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/gpu_env_test.log
# Verifies pb_gpu_venv sees the GPU and times a fp32 matmul (reference point for the training-cost table in DESIGN 6).
PB=/well/blagden/users/eic789/PerturbBridge
$PB/envs/pb_gpu_venv/bin/python - <<'PY'
import torch, time
print("torch", torch.__version__, "cuda available", torch.cuda.is_available(), "device", torch.cuda.get_device_name(0) if torch.cuda.is_available() else None)
x = torch.randn(8192, 8192, device="cuda"); torch.cuda.synchronize(); t = time.time()
for _ in range(10): y = x @ x
torch.cuda.synchronize(); dt = (time.time() - t) / 10; print(f"fp32 8192x8192 matmul {dt*1000:.1f} ms  = {2*8192**3/dt/1e12:.1f} TFLOPS")
xc = torch.randn(8192, 8192); t = time.time(); yc = xc @ xc; print(f"same on CPU ({torch.get_num_threads()} threads): {(time.time()-t)*1000:.0f} ms")
print("DONE gpuenvtest")
PY
