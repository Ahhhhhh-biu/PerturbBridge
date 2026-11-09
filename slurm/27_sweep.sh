#!/bin/bash
#SBATCH --job-name=pb_sweep
#SBATCH --partition=gpu_a100_80gb
#SBATCH --account=gpu_blagden.prj
#SBATCH --qos=gpu_bmrc_4hr
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --array=0-63
#SBATCH --output=/well/blagden/users/eic789/PerturbBridge/logs/sweep_%a.log
# Training-optimisation sweep (2026-11-07), tier A only: 8 configs x 2 cell lines x folds {0,1} x seeds {0,1} = 64 tasks.
# id = cfg*8 + cell*4 + fold*2 + seed (round 1 ids 0-63 = cfg 0-7; round 2 ids 64-127 = cfg 8-15). Output dirs: 13_model/runs/<cell>_f<fold>_A_sw<cfg>_s<seed>/
PB=/well/blagden/users/eic789/PerturbBridge
export TMPDIR=$PB/envs/tmp
id=$SLURM_ARRAY_TASK_ID; cfg=$((id / 8)); rest=$((id % 8)); cell=$((rest / 4)); fold=$(( (rest % 4) / 2 )); seed=$((rest % 2))
cells=(HepG2 K562); cl=${cells[$cell]}
case $cfg in
  0) opts="" ;;                                                              # first-round baseline
  1) opts="--sel weighted" ;;                                                # selection metric only
  2) opts="--sel weighted --sigw 10" ;;                                      # + significant-row loss weight
  3) opts="--sel weighted --sigw 10 --rpca 32 --rdrop 0.3" ;;               # + compressed, regularised RBP branch
  4) opts="--sel weighted --sigw 10 --cos --lr 3e-4 --epochs 20" ;;         # + cosine schedule, lower lr, longer
  5) opts="--sel weighted --sigw 10 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20" ;;   # all
  6) opts="--sel weighted --sigw 30 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20" ;;   # all, stronger loss weight
  7) opts="--sel weighted --sigw 10 --rpca 64 --rdrop 0.3 --cos --lr 3e-4 --epochs 20" ;;   # all, wider PCA
  # round 2 (ids 64-127): around config 6, larger significant-row weights, lr / weight decay / width / batch
  8) opts="--sel weighted --sigw 30 --cos --lr 3e-4 --epochs 20" ;;                                      # cfg 6 without PCA/dropout
  9) opts="--sel weighted --sigw 100 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20" ;;              # sigw 100
 10) opts="--sel weighted --sigw 300 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20" ;;              # sigw 300
 11) opts="--sel weighted --sigw 1000 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20" ;;             # ~significant rows only
 12) opts="--sel weighted --sigw 100 --rpca 32 --rdrop 0.3 --cos --lr 1e-4 --epochs 20" ;;              # sigw 100, lower lr
 13) opts="--sel weighted --sigw 100 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20 --wd 1e-2" ;;    # sigw 100, strong weight decay
 14) opts="--sel weighted --sigw 100 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20 --d 64" ;;       # sigw 100, narrower
 15) opts="--sel weighted --sigw 100 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20 --batch 2048" ;; # sigw 100, smaller batch
  # round 3 (ids 128-159): around cfg 15; calibration on null rows is the second criterion (cfg 11 inflated 16-22 % of null rows above 0.1)
 16) opts="--sel weighted --sigw 300 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20 --batch 2048" ;;            # sigw 300, batch 2048
 17) opts="--sel weighted --sigw 100 --rpca 32 --rdrop 0.3 --cos --lr 1e-4 --epochs 20 --batch 2048" ;;            # sigw 100, batch 2048, lr 1e-4
 18) opts="--sel weighted --sigw 1000 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20 --batch 2048" ;;           # sigw 1000, batch 2048
 19) opts="--sel weighted --sigw 100 --rpca 32 --rdrop 0.3 --cos --lr 3e-4 --epochs 20 --batch 1024" ;;            # sigw 100, batch 1024
esac
echo "task $id cfg $cfg $cl fold $fold seed $seed opts: $opts"
$PB/envs/pb_gpu_venv/bin/python $PB/13_model/01_train.py $cl $fold A --seed $seed --tag _sw${cfg}_s${seed} $opts && echo "DONE sweep $cl $fold cfg $cfg seed $seed $(date)"
