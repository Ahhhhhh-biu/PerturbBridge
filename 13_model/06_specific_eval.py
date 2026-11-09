"""RBP-specific evaluation: remove the event-common response and score what is left.
For each held-out fold, the event consensus ybar_e = mean of control-centred dPSI over TRAINING runs (runs not in the held-out fold).
Residual target  dy = y - ybar_e ; residual prediction dp = p - pbar_e where pbar_e = mean of the model's predictions for event e over the
held-out runs of that fold (what the model thinks is common to all RBPs). Scored on the same significant rows as sign_sig:
  sign_common   : sign(ybar_e) == sign(y)          -> how far event consensus alone gets (no model)
  sign_resid    : sign(dp) == sign(dy)             -> direction of the RBP-specific part
  sign_resid_big: same, rows with |dy| > 0.1       -> RBP-specific part that is large
  r_resid       : corr(dp, dy) on significant rows
Per run -> per RBP mean -> RBP mean per (version, cell, tier). Versions: v2 (_v2_s0-2), ctrl_none, ctrl_shuf.
Outputs 13_model/specific_eval_by_run.tsv, specific_eval.tsv (+ stdout)."""
import glob, re, numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; out = []
for cl in ("HepG2", "K562"):
    z = np.load(f"{PB}/13_model/tensors_{cl}.npz", allow_pickle=True); e_all, r_all, y_all, sig_all = z["row_e"], z["row_r"], z["row_y"], z["row_sig"].astype(bool)
    rf = z["run_fold"][r_all]; ne = int(e_all.max()) + 1; targets = z["targets"]; run_target = z["run_target"]
    ybar = {}
    for fold in range(5):
        tr = rf != fold; s = np.bincount(e_all[tr], weights=y_all[tr], minlength=ne); c = np.bincount(e_all[tr], minlength=ne); ybar[fold] = s / np.maximum(c, 1)
    for f in sorted(glob.glob(f"{PB}/13_model/runs/{cl}_f[0-4]_[ABC]_v2_s[0-9]/pred_heldout.npz") + glob.glob(f"{PB}/13_model/runs/{cl}_f[0-4]_A_ctrl_*_s[0-9]/pred_heldout.npz")):
        dn = f.split("/")[-2]; m = re.match(rf"{cl}_f(\d)_([ABC])_(v2|ctrl_none|ctrl_shuf)_s(\d)$", dn)
        if not m: continue
        fold, tier, ver, seed = int(m.group(1)), m.group(2), m.group(3), int(m.group(4)); d = np.load(f); idx, p = d["idx"], d["pred"]
        e, r, y, sig = e_all[idx], r_all[idx], y_all[idx], sig_all[idx]
        ps = np.bincount(e, weights=p, minlength=ne); pc = np.bincount(e, minlength=ne); pbar = ps / np.maximum(pc, 1)
        yb = ybar[fold][e]; dy = y - yb; dp = p - pbar[e]
        df = pd.DataFrame({"r": r, "y": y, "yb": yb, "dy": dy, "dp": dp, "p": p, "sig": sig})
        for run, g in df.groupby("r"):
            gs = g[g.sig]
            if len(gs) < 20: continue
            big = gs[gs.dy.abs() > 0.1]
            out.append({"version": ver, "cell_line": cl, "tier": tier, "fold": fold, "seed": seed, "run": int(run), "target": targets[run_target[run]], "n_sig": len(gs),
                        "sign_sig": (np.sign(gs.p) == np.sign(gs.y)).mean(), "sign_common": (np.sign(gs.yb) == np.sign(gs.y)).mean(),
                        "sign_resid": (np.sign(gs.dp) == np.sign(gs.dy)).mean(), "n_resid_big": len(big), "sign_resid_big": (np.sign(big.dp) == np.sign(big.dy)).mean() if len(big) >= 10 else np.nan,
                        "r_resid": np.corrcoef(gs.dp, gs.dy)[0, 1] if gs.dp.std() > 0 else np.nan, "r_sig": np.corrcoef(gs.p, gs.y)[0, 1] if gs.p.std() > 0 else np.nan,
                        "frac_resid_big": (gs.dy.abs() > 0.1).mean()})
d = pd.DataFrame(out); d.to_csv(f"{PB}/13_model/specific_eval_by_run.tsv", sep="\t", index=False)
met = ["sign_sig", "sign_common", "sign_resid", "sign_resid_big", "r_sig", "r_resid", "frac_resid_big"]
byr = d.groupby(["version", "cell_line", "tier", "target"])[met].mean().reset_index()
g = byr.groupby(["version", "cell_line", "tier"]).agg(n_rbp=("target", "size"), **{m: (m, "mean") for m in met}).reset_index(); g.to_csv(f"{PB}/13_model/specific_eval.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(g.round(3).to_string(index=False))
# paired v2 (tier A) minus ctrl_none on same RBPs, bootstrap CI
rng = np.random.default_rng(0)
for cl in ("HepG2", "K562"):
    a = byr[(byr.version == "v2") & (byr.cell_line == cl) & (byr.tier == "A")].set_index("target"); b = byr[(byr.version == "ctrl_none") & (byr.cell_line == cl)].set_index("target")
    for m in ("sign_sig", "sign_resid", "r_sig", "r_resid"):
        dd = (a[m] - b[m].reindex(a.index)).dropna().values; bs = [rng.choice(dd, len(dd)).mean() for _ in range(2000)]
        print(f"{cl} v2(A) - event-only, {m}: {dd.mean():+.4f} [{np.percentile(bs, 2.5):+.4f},{np.percentile(bs, 97.5):+.4f}] n={len(dd)}")
