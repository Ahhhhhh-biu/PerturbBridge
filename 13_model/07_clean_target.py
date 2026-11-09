"""Prototype of a control-independent, replicate-consistent target set ("clean set").
Significant events in ENCODE are selected by KD vs ITS OWN control; a control sample that happens to sit far from the event's typical PSI
produces a 'significant' dPSI in the opposite direction (regression to the mean). A rule with no RBP information, sign(-(psi_ctrl - typical PSI)),
gets 0.83/0.79 direction accuracy on those rows (logs/ctrl_artefact_check.txt). This script defines events by the KD sample alone:
  consensus_e = leave-one-out mean of KD PSI over the other runs sharing the same control  (= psi_kd - dpsi_cc, by construction)
  d1 = kd_rep1 - consensus, d2 = kd_rep2 - consensus
  clean = sign(d1) == sign(d2) & min(|d1|, |d2|) > 0.05 & |dpsi_cc| > 0.1   (strict: min > 0.1)
and scores on the clean set: the control-deviation rule, model v2 (tier A, seed 0), the event-only control and the shuffled-RBP control.
Row order of the tensors is checked against the parquet. Outputs 13_model/clean_target_summary.tsv, clean_target_by_rbp.tsv (+ stdout)."""
import glob, numpy as np, pandas as pd, pyarrow.parquet as pq
PB = "/well/blagden/users/eic789/PerturbBridge"; U = f"{PB}/10_event_universe/encode_full"; rng = np.random.default_rng(0)
def acc(p, y): return (np.sign(p) == np.sign(y)).mean()
S, B = [], []
for cl in ("HepG2", "K562"):
    t = pq.read_table(f"{U}/train_table_{cl}.parquet", columns=["event_id", "run", "psi_kd", "psi_ctrl", "dpsi", "dpsi_cc", "fdr", "psi_kd_reps"]).to_pandas()
    t = t[np.isfinite(t.dpsi_cc)].reset_index(drop=True)
    z = np.load(f"{PB}/13_model/tensors_{cl}.npz", allow_pickle=True); y = z["row_y"]; assert len(t) == len(y) and np.allclose(t.dpsi_cc.values, y, atol=1e-5), "row order mismatch"
    e, r, psi, sig = z["row_e"], z["row_r"], z["row_psi"], z["row_sig"].astype(bool); ne = e.max() + 1
    pm = np.bincount(e, weights=psi, minlength=ne) / np.bincount(e, minlength=ne); dev = psi - pm[e]
    sp = t.psi_kd_reps.str.split(",", n=2, expand=True); k1 = pd.to_numeric(sp[0], errors="coerce").values; k2 = pd.to_numeric(sp[1], errors="coerce").values
    cons = t.psi_kd.values - y; d1 = k1 - cons; d2 = k2 - cons; both = np.isfinite(d1 * d2)
    same = both & (np.sign(d1) == np.sign(d2)); mn = np.minimum(np.abs(d1), np.abs(d2))
    sets = {"sig_FDR": sig, "clean_0.05": same & (mn > 0.05) & (np.abs(y) > 0.1), "clean_0.1": same & (mn > 0.1) & (np.abs(y) > 0.1), "clean_0.05_notsig": same & (mn > 0.05) & (np.abs(y) > 0.1) & ~sig,
            "cc>0.1_inconsistent": both & (np.abs(y) > 0.1) & ~(same & (mn > 0.05))}
    P = {}
    for ver, pat in (("v2_A", f"{cl}_f[0-4]_A_v2_s0"), ("event_only", f"{cl}_f[0-4]_A_ctrl_none_s0"), ("shuffled", f"{cl}_f[0-4]_A_ctrl_shuf_s0")):
        pr = np.full(len(y), np.nan, dtype=np.float32)
        for f in glob.glob(f"{PB}/13_model/runs/{pat}/pred_heldout.npz"): d = np.load(f); pr[d["idx"]] = d["pred"]
        P[ver] = pr
    ok = np.isfinite(P["v2_A"]) & np.isfinite(P["event_only"]) & np.isfinite(P["shuffled"])
    print(f"== {cl}: rows {len(y):,} | with 2 KD reps {both.mean():.3f} ==")
    print(f"{'row set':22s} {'n':>10s} {'frac':>6s} {'overlap sig':>11s} {'-ctrl_dev rule':>15s} {'v2_A':>7s} {'event_only':>11s} {'shuffled':>9s} {'mean|dev|':>9s}")
    for nm, m in sets.items():
        m = m & ok; row = {"cell_line": cl, "set": nm, "n": int(m.sum()), "frac": m.mean(), "overlap_sig": (m & sig).sum() / max(m.sum(), 1), "rule": acc(-dev[m], y[m]), "v2_A": acc(P["v2_A"][m], y[m]), "event_only": acc(P["event_only"][m], y[m]), "shuffled": acc(P["shuffled"][m], y[m]), "mean_abs_dev": np.abs(dev[m]).mean()}
        S.append(row); print(f"{nm:22s} {row['n']:10,d} {row['frac']:6.3f} {row['overlap_sig']:11.3f} {row['rule']:15.3f} {row['v2_A']:7.3f} {row['event_only']:11.3f} {row['shuffled']:9.3f} {row['mean_abs_dev']:9.3f}")
    # per RBP on clean_0.05: RBP-mean accuracy and paired differences
    m = sets["clean_0.05"] & ok; df = pd.DataFrame({"r": r[m], "y": y[m], **{k: v[m] for k, v in P.items()}, "rule": -dev[m]})
    tg = z["targets"][z["run_target"]]
    g = df.groupby("r").apply(lambda d: pd.Series({"n": len(d), **{k: acc(d[k], d.y) for k in ("v2_A", "event_only", "shuffled", "rule")}})); g["target"] = tg[g.index]; g = g[g.n >= 20]
    br = g.groupby("target").mean(numeric_only=True); br["cell_line"] = cl; B.append(br.reset_index())
    print(f"per-RBP mean on clean_0.05 (runs with >=20 clean rows, {len(br)} RBPs): v2_A {br.v2_A.mean():.3f} | event_only {br.event_only.mean():.3f} | shuffled {br.shuffled.mean():.3f} | rule {br.rule.mean():.3f} | median clean rows per run {g.n.median():.0f}")
    for a_, b_ in (("v2_A", "event_only"), ("v2_A", "shuffled"), ("v2_A", "rule")):
        dd = (br[a_] - br[b_]).values; bs = [rng.choice(dd, len(dd)).mean() for _ in range(2000)]; print(f"   {a_} - {b_}: {dd.mean():+.4f} [{np.percentile(bs, 2.5):+.4f},{np.percentile(bs, 97.5):+.4f}]")
    for nm in ("SF3B1", "SF3A1", "U2AF2", "PTBP1", "HNRNPC", "QKI", "SRSF1", "LARP4", "RBFOX2"):
        if nm in br.index: print(f"   {nm}: v2 {br.loc[nm, 'v2_A']:.3f} | event_only {br.loc[nm, 'event_only']:.3f} | shuffled {br.loc[nm, 'shuffled']:.3f} | rule {br.loc[nm, 'rule']:.3f}")
    del t
pd.DataFrame(S).to_csv(f"{PB}/13_model/clean_target_summary.tsv", sep="\t", index=False); pd.concat(B).to_csv(f"{PB}/13_model/clean_target_by_rbp.tsv", sep="\t", index=False)
