"""Summarise the training-optimisation sweep (slurm/27_sweep.sh): runs/<cell>_f<fold>_A_sw<cfg>_s<seed>/.
For each config, paired difference vs config 0 on the same (cell, fold, seed, RBP): RBP-mean and n_sig-weighted sign_sig and r_sig,
bootstrap CI over RBPs. Also the epoch-to-epoch instability of the inner-validation score (std of per-epoch val r_sig / its mean range)
from log.txt. Output: 13_model/sweep_summary.tsv (+ stdout)."""
import glob, re, warnings, numpy as np, pandas as pd; warnings.simplefilter("ignore")
PB = "/well/blagden/users/eic789/PerturbBridge"; rng = np.random.default_rng(0)
CFG = {16: "sigw300 b2048", 17: "sigw100 b2048 lr1e-4", 18: "sigw1000 b2048", 19: "sigw100 b1024", 0: "baseline", 1: "sel weighted", 2: "+sigw 10", 3: "+rpca32 rdrop0.3", 4: "+cos lr3e-4 20ep", 5: "all", 6: "all sigw30", 7: "all rpca64", 8: "sigw30 cos noPCA", 9: "all sigw100", 10: "all sigw300", 11: "all sigw1000", 12: "sigw100 lr1e-4", 13: "sigw100 wd1e-2", 14: "sigw100 d64", 15: "sigw100 batch2048"}
R, V = [], []
for f in sorted(glob.glob(f"{PB}/13_model/runs/*_f[0-9]_A_sw[0-9]*_s[0-9]/per_run.tsv")):
    dn = f.split("/")[-2]; m = re.match(r"(\w+)_f(\d)_A_sw(\d+)_s(\d)", dn); d = pd.read_csv(f, sep="\t")
    d["cell_line"], d["fold"], d["cfg"], d["seed"] = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4)); R.append(d)
    ep = [float(x) for x in re.findall(r"val r_sig\(by RBP\) (-?[\d.]+)", open(f.replace("per_run.tsv", "log.txt")).read())]
    best = re.findall(r"best epoch (\d+)", open(f.replace("per_run.tsv", "log.txt")).read())
    V.append(dict(cell_line=m.group(1), fold=int(m.group(2)), cfg=int(m.group(3)), seed=int(m.group(4)), n_ep=len(ep), val_mean=np.mean(ep), val_std=np.std(ep), val_last3=np.mean(ep[-3:]), best_ep=int(best[-1]) if best else -1))
R = pd.concat(R, ignore_index=True); V = pd.DataFrame(V)
print("tasks found:", R.groupby(["cfg"]).apply(lambda g: g.groupby(["cell_line", "fold", "seed"]).ngroups).to_dict())
by = R.groupby(["cell_line", "fold", "seed", "cfg", "target"]).agg(r_sig=("r_sig", "mean"), sign_sig=("sign_sig", "mean"), n_sig=("n_sig", "sum")).reset_index()
def paired(d, w):
    bs = [np.average(d[i := rng.integers(0, len(d), len(d))], weights=w[i]) for _ in range(2000)]
    return f"{np.average(d, weights=w):+.4f} [{np.percentile(bs, 2.5):+.4f},{np.percentile(bs, 97.5):+.4f}]"
out = []
for cl in ("HepG2", "K562"):
    b = by[by.cell_line == cl]; base = b[b.cfg == 0].set_index(["fold", "seed", "target"])
    for cfg in sorted(b.cfg.unique()):
        c = b[b.cfg == cfg].set_index(["fold", "seed", "target"]); j = c.join(base, rsuffix="_0", how="inner").dropna(subset=["sign_sig", "sign_sig_0"])
        w = j.n_sig.values.astype(float); one = np.ones(len(j)); jr = j.dropna(subset=["r_sig", "r_sig_0"]); wr = jr.n_sig.values.astype(float)
        vv = V[(V.cell_line == cl) & (V.cfg == cfg)]
        row = dict(cell_line=cl, cfg=cfg, name=CFG.get(cfg, ""), n_pairs=len(j), sign_sig=f"{j.sign_sig.mean():.4f}", sign_sig_w=f"{np.average(j.sign_sig, weights=w):.4f}", r_sig_w=f"{np.average(jr.r_sig, weights=wr):.4f}",
                   d_sign=paired((j.sign_sig - j.sign_sig_0).values, one) if cfg else "", d_sign_w=paired((j.sign_sig - j.sign_sig_0).values, w) if cfg else "",
                   d_rsig_w=paired((jr.r_sig - jr.r_sig_0).values, wr) if cfg else "", val_std=f"{vv.val_std.mean():.4f}", best_ep=",".join(map(str, vv.best_ep.tolist())))
        out.append(row)
S = pd.DataFrame(out); S.to_csv(f"{PB}/13_model/sweep_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 20); print(S.to_string(index=False))
