"""Measurement ceilings for direction accuracy, computed without selection bias, to compare with the model's held-out sign_sig.
(a) within-run split-half: D1 = kd_rep1 - ctrl_rep1, D2 = kd_rep2 - ctrl_rep2; select events with |D1| > thr on half 1 only, score
    sign agreement of D2 on them (and the symmetric direction); per run, then per RBP. Single replicate vs single replicate, so this
    is a LOWER bound for the two-replicate consensus the model is scored against.
(b) cross-run, same RBP, same cell line (shRNA vs CRISPR or two shRNA runs): select significant events in run A (FDR<0.05 & |dPSI|>0.1,
    exactly the model's test set), score sign agreement of control-centred dPSI in run B. An independent experiment predicting the
    model's target: the fairest ceiling, but only ~25 pairs and mostly across technologies.
Also: per-RBP within-run reproducibility vs model v2 tier-A sign_sig (results_by_rbp_v2.tsv) -> scatter.
Outputs: 13_model/ceiling_within_run.tsv, ceiling_within_run_by_rbp.tsv, ceiling_cross_run.tsv, fig_ceiling_vs_model.png (300 dpi); stdout summary."""
import numpy as np, pandas as pd, pyarrow.parquet as pq, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PB = "/well/blagden/users/eic789/PerturbBridge"; U = f"{PB}/10_event_universe/encode_full"
def first2(col):
    sp = col.str.split(",", n=2, expand=True); return pd.to_numeric(sp[0], errors="coerce").values, pd.to_numeric(sp[1], errors="coerce").values
W, X = [], []
for cl in ("HepG2", "K562"):
    t = pq.read_table(f"{U}/train_table_{cl}.parquet", columns=["event_id", "run", "target", "assay", "dpsi", "dpsi_cc", "fdr", "psi_kd_reps", "psi_ctrl_reps"]).to_pandas()
    k1, k2 = first2(t.psi_kd_reps); c1, c2 = first2(t.psi_ctrl_reps); t["D1"] = k1 - c1; t["D2"] = k2 - c2; t = t.drop(columns=["psi_kd_reps", "psi_ctrl_reps"])
    t["sig"] = (t.fdr < 0.05) & (t.dpsi.abs() > 0.1)
    print(cl, "rows", len(t), "| rows with both halves", int(np.isfinite(t.D1 * t.D2).sum()), "| runs with 2 reps", t[np.isfinite(t.D1 * t.D2)].run.nunique(), "/", t.run.nunique())
    for run, g in t.groupby("run"):
        ok = np.isfinite(g.D1 * g.D2); g = g[ok]
        if len(g) < 1000: continue
        row = {"cell_line": cl, "run": run, "target": g.target.iloc[0], "assay": g.assay.iloc[0], "n": len(g), "n_sig": int(g.sig.sum())}
        for thr in (0.1, 0.2):
            s1 = g.D1.abs() > thr; s2 = g.D2.abs() > thr
            a12 = (np.sign(g.D2[s1]) == np.sign(g.D1[s1])).mean() if s1.sum() >= 20 else np.nan
            a21 = (np.sign(g.D1[s2]) == np.sign(g.D2[s2])).mean() if s2.sum() >= 20 else np.nan
            row[f"agree_{thr}"] = np.nanmean([a12, a21]); row[f"n_sel_{thr}"] = int((s1.sum() + s2.sum()) / 2)
        # biased reference: direction agreement of the two halves on the FULL-run significant events (selected with both halves)
        sg = g[g.sig]; row["agree_on_fullsig_biased"] = (np.sign(sg.D1) == np.sign(sg.D2)).mean() if len(sg) >= 20 else np.nan
        row["r_all_halves"] = np.corrcoef(g.D1, g.D2)[0, 1]
        W.append(row)
    # (b) cross-run pairs
    runs = t[["run", "target", "assay"]].drop_duplicates(); multi = runs.groupby("target").run.nunique(); multi = multi[multi >= 2].index
    for tg in multi:
        rr = runs[runs.target == tg].run.tolist()
        for ra in rr:
            for rb in rr:
                if ra == rb: continue
                A = t[(t.run == ra) & t.sig][["event_id", "dpsi", "dpsi_cc"]]; B = t[t.run == rb][["event_id", "dpsi_cc", "sig"]]
                m = A.merge(B, on="event_id", suffixes=("_a", "_b"))
                if len(m) < 20: continue
                X.append({"cell_line": cl, "target": tg, "run_select": ra, "assay_select": runs.set_index("run").assay[ra], "run_test": rb, "assay_test": runs.set_index("run").assay[rb],
                          "n_sig_a": len(A), "n_overlap": len(m), "sign_agree_cc": (np.sign(m.dpsi_cc_a) == np.sign(m.dpsi_cc_b)).mean(), "sign_agree_raw": (np.sign(m.dpsi) == np.sign(m.dpsi_cc_b)).mean(),
                          "r_sig_cc": np.corrcoef(m.dpsi_cc_a, m.dpsi_cc_b)[0, 1], "frac_also_sig_b": m.sig.mean()})
    del t
W = pd.DataFrame(W); W.to_csv(f"{PB}/13_model/ceiling_within_run.tsv", sep="\t", index=False)
X = pd.DataFrame(X); X.to_csv(f"{PB}/13_model/ceiling_cross_run.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print("\n== (a) within-run split-half, per cell line: mean over runs (median) ==")
for cl, g in W.groupby("cell_line"):
    print(cl, f"runs {len(g)} | agree_0.1 {g['agree_0.1'].mean():.3f} ({g['agree_0.1'].median():.3f}) | agree_0.2 {g['agree_0.2'].mean():.3f} ({g['agree_0.2'].median():.3f}) | biased-on-fullsig {g.agree_on_fullsig_biased.mean():.3f} | r_all halves {g.r_all_halves.median():.3f} | n_sel_0.1 median {g['n_sel_0.1'].median():.0f}")
    print("   by assay:", g.groupby("assay")["agree_0.1"].agg(["mean", "count"]).round(3).to_dict())
    big = g[g.n_sig >= 1000]; print(f"   runs with >=1000 sig events: {len(big)} agree_0.1 {big['agree_0.1'].mean():.3f}")
by = W.groupby(["cell_line", "target"]).agg(agree_0_1=("agree_0.1", "mean"), agree_0_2=("agree_0.2", "mean"), n_sig=("n_sig", "sum"), n_runs=("run", "size")).reset_index()
v2 = pd.read_csv(f"{PB}/13_model/results_by_rbp_v2.tsv", sep="\t"); v2 = v2[v2.method == "model_A"][["cell_line", "target", "sign_sig", "r_sig"]]
by = by.merge(v2, on=["cell_line", "target"], how="left"); by.to_csv(f"{PB}/13_model/ceiling_within_run_by_rbp.tsv", sep="\t", index=False)
print("\n== per-RBP: within-run reproducibility vs model v2 (tier A) ==")
fig, ax = plt.subplots(1, 2, figsize=(8, 4), dpi=300)
for i, cl in enumerate(("HepG2", "K562")):
    g = by[(by.cell_line == cl) & by.sign_sig.notna() & by.agree_0_1.notna()]
    r = np.corrcoef(g.agree_0_1, g.sign_sig)[0, 1]; above = (g.sign_sig > g.agree_0_1).mean()
    print(f"{cl}: n_rbp {len(g)} | replicate agreement mean {g.agree_0_1.mean():.3f} | model sign_sig mean {g.sign_sig.mean():.3f} | r(agreement, model) {r:.3f} | RBPs where model > replicate agreement {above:.2f}")
    for nm in ("SF3B1", "U2AF2", "PTBP1", "HNRNPC", "QKI", "RBFOX2", "SRSF1", "LARP4"):
        s = g[g.target == nm]
        if len(s): print(f"   {nm}: replicate {s.agree_0_1.iloc[0]:.3f} | model {s.sign_sig.iloc[0]:.3f} | n_sig {int(s.n_sig.iloc[0])}")
    ax[i].scatter(g.agree_0_1, g.sign_sig, s=np.clip(g.n_sig / 50, 4, 80), alpha=0.5, lw=0); ax[i].plot([0.4, 1], [0.4, 1], "k--", lw=0.6)
    ax[i].axhline(0.58, color="grey", lw=0.6, ls=":"); ax[i].set_xlabel("within-run replicate direction agreement (|ΔPSI|>0.1 in one half)"); ax[i].set_ylabel("model v2 held-out sign_sig (tier A)"); ax[i].set_title(f"{cl}  r={r:.2f}  n={len(g)}", fontsize=9)
plt.tight_layout(); plt.savefig(f"{PB}/13_model/fig_ceiling_vs_model.png"); print("fig", f"{PB}/13_model/fig_ceiling_vs_model.png")
print("\n== (b) cross-run same-RBP pairs: select sig in run A, test control-centred sign in run B ==")
print(X.round(3).to_string(index=False))
for cl, g in X.groupby("cell_line"): print(cl, f"pairs {len(g)} | sign_agree_cc mean {g.sign_agree_cc.mean():.3f} median {g.sign_agree_cc.median():.3f} | cross-technology only: {g[g.assay_select != g.assay_test].sign_agree_cc.mean():.3f} (n={int((g.assay_select != g.assay_test).sum())}) | r_sig_cc mean {g.r_sig_cc.mean():.3f}")
# model v2 on the same test runs (per_run files, tier A, mean over seeds)
import glob
pr = pd.concat([pd.read_csv(f, sep="\t").assign(cl=f.split("/")[-2].split("_")[0]) for f in glob.glob(f"{PB}/13_model/runs/*_f[0-4]_A_v2_s[0-2]/per_run.tsv")])
mr = pr.groupby(["cl", "run"]).sign_sig.mean().rename("model_sign_sig").reset_index().rename(columns={"cl": "cell_line", "run": "run_test"})
X2 = X.merge(mr, on=["cell_line", "run_test"], how="left")
print("\nmodel v2 (tier A) on the same test runs vs replicate-experiment ceiling:")
print(X2[["cell_line", "target", "assay_select", "assay_test", "n_overlap", "sign_agree_cc", "model_sign_sig"]].round(3).to_string(index=False))
print(f"mean over pairs: replicate experiment {X2.sign_agree_cc.mean():.3f} | model {X2.model_sign_sig.mean():.3f}")
