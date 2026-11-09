"""Calibration of held-out predictions: for each run dir (13_model/runs/<cell>_f<k>_<tier>[<tag>]_s<seed>/pred_heldout.npz) compute, over
held-out rows, the fraction of non-significant rows predicted |pred|>0.1 (false alarms), the fraction of significant rows predicted
|pred|>0.1 (hits), mean |pred| on each, and MSE; then average over folds and seeds per (version, cell, tier).
Usage: python 04_calibration.py  -> 13_model/calibration.tsv (+ stdout). Versions: v1 = runs without tag, v2 = tag _v2."""
import glob, re, numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; out = []
T = {cl: np.load(f"{PB}/13_model/tensors_{cl}.npz") for cl in ("HepG2", "K562")}
for f in sorted(glob.glob(f"{PB}/13_model/runs/*_f[0-4]_[ABC]/pred_heldout.npz") + glob.glob(f"{PB}/13_model/runs/*_f[0-4]_[ABC]_s[0-9]/pred_heldout.npz") + glob.glob(f"{PB}/13_model/runs/*_f[0-4]_[ABC]_v2_s[0-9]/pred_heldout.npz")):
    dn = f.split("/")[-2]; m = re.match(r"(\w+?)_f(\d)_([ABC])(_v2)?(?:_s(\d))?$", dn)
    if not m: continue
    cl, fold, tier, tag, seed = m.groups(); z = np.load(f); idx, p = z["idx"], z["pred"]
    y = T[cl]["row_y"][idx]; sig = T[cl]["row_sig"][idx].astype(bool); ap = np.abs(p)
    out.append({"version": "v2" if tag else "v1", "cell_line": cl, "tier": tier, "fold": int(fold), "seed": int(seed or 0), "n_rows": len(p), "frac_sig": sig.mean(),
                "null_gt0.1": (ap[~sig] > 0.1).mean(), "sig_gt0.1": (ap[sig] > 0.1).mean(), "mean_abs_pred_null": ap[~sig].mean(), "mean_abs_pred_sig": ap[sig].mean(),
                "mean_abs_y_sig": np.abs(y[sig]).mean(), "mse_all": ((p - y) ** 2).mean(), "mse_sig": ((p[sig] - y[sig]) ** 2).mean(), "mse_null": ((p[~sig] - y[~sig]) ** 2).mean()})
d = pd.DataFrame(out); d.to_csv(f"{PB}/13_model/calibration_by_run.tsv", sep="\t", index=False)
g = d.groupby(["version", "cell_line", "tier"]).agg(n=("fold", "size"), **{c: (c, "mean") for c in d.columns[6:]}).reset_index(); g.to_csv(f"{PB}/13_model/calibration.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(g.round(4).to_string(index=False))
