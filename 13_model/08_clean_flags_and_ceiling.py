"""(1) Write the control-independent replicate-consistent 'clean' flag per training row (13_model/clean_flags_<cell>.npz: row_clean, row_cons =
leave-one-out consensus KD PSI, aligned to tensors_<cell>.npz rows; definition as in 07_clean_target.py, min |d| > 0.05, |dpsi_cc| > 0.1).
(2) Ceilings on the clean set: (a) cross-run same-RBP pairs: clean events of run A -> sign of dpsi_cc in run B; (b) zero point: for random pairs
of DIFFERENT RBPs (same cell line), direction agreement on events clean in both runs (the event-common response). Outputs
13_model/clean_ceiling_cross_run.tsv, clean_zero_point.tsv (+ stdout)."""
import numpy as np, pandas as pd, pyarrow.parquet as pq
PB = "/well/blagden/users/eic789/PerturbBridge"; U = f"{PB}/10_event_universe/encode_full"; rng = np.random.default_rng(0); X, Z = [], []
for cl in ("HepG2", "K562"):
    t = pq.read_table(f"{U}/train_table_{cl}.parquet", columns=["event_id", "run", "target", "assay", "psi_kd", "dpsi_cc", "psi_kd_reps"]).to_pandas()
    t = t[np.isfinite(t.dpsi_cc)].reset_index(drop=True); z = np.load(f"{PB}/13_model/tensors_{cl}.npz", allow_pickle=True); y = z["row_y"]
    assert len(t) == len(y) and np.allclose(t.dpsi_cc.values, y, atol=1e-5)
    sp = t.psi_kd_reps.str.split(",", n=2, expand=True); k1 = pd.to_numeric(sp[0], errors="coerce").values; k2 = pd.to_numeric(sp[1], errors="coerce").values
    cons = t.psi_kd.values - y; d1 = k1 - cons; d2 = k2 - cons
    clean = np.isfinite(d1 * d2) & (np.sign(d1) == np.sign(d2)) & (np.minimum(np.abs(d1), np.abs(d2)) > 0.05) & (np.abs(y) > 0.1)
    np.savez_compressed(f"{PB}/13_model/clean_flags_{cl}.npz", row_clean=clean, row_cons=cons.astype(np.float32)); print(cl, "clean rows", clean.sum(), f"{clean.mean():.4f}", "| per run median", pd.Series(clean).groupby(t.run.values).sum().median())
    t["clean"] = clean; t = t[["event_id", "run", "target", "assay", "dpsi_cc", "clean"]]
    runs = t[["run", "target", "assay"]].drop_duplicates().set_index("run"); multi = runs.groupby("target").size(); multi = multi[multi >= 2].index
    byrun = {run: g.set_index("event_id") for run, g in t[t.run.isin(runs[runs.target.isin(multi)].index)].groupby("run")}
    for tg in multi:
        rr = runs[runs.target == tg].index.tolist()
        for ra in rr:
            for rb in rr:
                if ra == rb: continue
                A = byrun[ra]; A = A[A.clean]; B = byrun[rb]; m = A.join(B[["dpsi_cc", "clean"]], rsuffix="_b", how="inner")
                if len(m) < 20: continue
                X.append({"cell_line": cl, "target": tg, "run_select": ra, "assay_select": runs.assay[ra], "run_test": rb, "assay_test": runs.assay[rb], "n_clean_a": len(A), "n_overlap": len(m),
                          "sign_agree": (np.sign(m.dpsi_cc) == np.sign(m.dpsi_cc_b)).mean(), "sign_agree_both_clean": (np.sign(m.dpsi_cc[m.clean_b]) == np.sign(m.dpsi_cc_b[m.clean_b])).mean() if m.clean_b.sum() >= 20 else np.nan, "frac_also_clean_b": m.clean_b.mean()})
    # zero point: 400 random pairs of runs with different targets
    allruns = runs.index.tolist(); tc = t[t.clean].groupby("run")
    groups = {run: g.set_index("event_id").dpsi_cc for run, g in tc}
    for _ in range(400):
        ra, rb = rng.choice(allruns, 2, replace=False)
        if runs.target[ra] == runs.target[rb] or ra not in groups or rb not in groups: continue
        m = groups[ra].to_frame("a").join(groups[rb].to_frame("b"), how="inner")
        if len(m) < 20: continue
        Z.append({"cell_line": cl, "run_a": ra, "run_b": rb, "n_both_clean": len(m), "sign_agree": (np.sign(m.a) == np.sign(m.b)).mean()})
    del t
X = pd.DataFrame(X); X.to_csv(f"{PB}/13_model/clean_ceiling_cross_run.tsv", sep="\t", index=False); Z = pd.DataFrame(Z); Z.to_csv(f"{PB}/13_model/clean_zero_point.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print("\n== clean-set ceiling: same RBP, independent run ==")
for cl, g in X.groupby("cell_line"): print(cl, f"pairs {len(g)} | sign_agree mean {g.sign_agree.mean():.3f} median {g.sign_agree.median():.3f} | on events clean in both runs {g.sign_agree_both_clean.mean():.3f} | frac also clean in B {g.frac_also_clean_b.mean():.3f} | cross-technology {g[g.assay_select != g.assay_test].sign_agree.mean():.3f}")
print(X.round(3).to_string(index=False))
print("\n== clean-set zero point: different RBPs, events clean in both ==")
for cl, g in Z.groupby("cell_line"): print(cl, f"pairs {len(g)} | sign_agree mean {g.sign_agree.mean():.3f} median {g.sign_agree.median():.3f} | n_both_clean median {g.n_both_clean.median():.0f}")
