"""H2 measurement ceilings on the junction-level training tables (design 4.3 positive control).
(1) split-half within run: dPSI from replicate 1 (KD rep1 - ctrl rep1) vs replicate 2, per run.
(2) cross-assay: same RBP, same cell line, shRNA vs CRISPR run, over events covered in both.
(3) null: dPSI correlation between runs of different RBPs in the same cell line, separately for pairs sharing / not sharing
    the control experiment (shared controls inflate the correlation).
Usage: python 01_ceiling.py  -> 12_baselines/ceiling_split_half.tsv, ceiling_cross_assay.tsv, ceiling_null_pairs.tsv (+ stdout)"""
import numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; U = f"{PB}/10_event_universe/encode_full"; O = f"{PB}/12_baselines"; rng = np.random.default_rng(1)
tiers = pd.read_csv(f"{PB}/11_splits/run_tiers.tsv", sep="\t")
def r(a, b):
    m = np.isfinite(a) & np.isfinite(b); return (np.corrcoef(a[m], b[m])[0, 1] if m.sum() > 2 else np.nan), int(m.sum())
sh, ca, nl = [], [], []
for cl in ["HepG2", "K562"]:
    t = pd.read_parquet(f"{U}/train_table_{cl}.parquet", columns=["event_id", "run", "target", "assay", "dpsi", "fdr", "psi_kd_reps", "psi_ctrl_reps"]).merge(tiers[["run", "control"]], on="run")
    k = t.psi_kd_reps.str.split(",", expand=True).iloc[:, :2].apply(pd.to_numeric, errors="coerce"); c = t.psi_ctrl_reps.str.split(",", expand=True).iloc[:, :2].apply(pd.to_numeric, errors="coerce")
    t["d1"] = k[0] - c[0]; t["d2"] = k[1] - c[1]; t["sig"] = (t.fdr < 0.05) & (t.dpsi.abs() > 0.1)
    for run, g in t.groupby("run"):
        ra, n = r(g.d1.values, g.d2.values); rs, ns = r(g.d1[g.sig].values, g.d2[g.sig].values)
        both = g[(g.d1.abs() > 0.1) & (g.d2.abs() > 0.1)]; sgn = (np.sign(both.d1) == np.sign(both.d2)).mean() if len(both) else np.nan
        sh.append((cl, run, g.target.iloc[0], g.assay.iloc[0], n, ra, ns, rs, len(both), sgn))
    M = t.pivot(index="event_id", columns="run", values="dpsi"); S = t.pivot(index="event_id", columns="run", values="sig").fillna(False)
    info = tiers.set_index("run").loc[M.columns]
    for tg, g in info.groupby("target"):
        if g.assay.nunique() == 2:
            a = g[g.assay.str.startswith("shRNA")].index[0]; b = g[g.assay.str.startswith("CRISPR")].index[0]
            ra, n = r(M[a].values, M[b].values); e = S[a] | S[b]; rs, ns = r(M[a][e].values, M[b][e].values)
            bb = S[a] & S[b]; sgn = (np.sign(M[a][bb]) == np.sign(M[b][bb])).mean() if bb.sum() else np.nan
            ca.append((cl, tg, a, b, n, ra, ns, rs, int(bb.sum()), sgn, g.loc[a, "control"] == g.loc[b, "control"]))
    runs = list(M.columns)
    for _ in range(3000):
        a, b = rng.choice(runs, 2, replace=False)
        if info.loc[a, "target"] == info.loc[b, "target"]: continue
        ra, n = r(M[a].values, M[b].values); e = S[a] | S[b]; rs, ns = r(M[a][e].values, M[b][e].values)
        nl.append((cl, a, b, info.loc[a, "control"] == info.loc[b, "control"], info.loc[a, "assay"] == info.loc[b, "assay"], n, ra, ns, rs))
sh = pd.DataFrame(sh, columns=["cell_line", "run", "target", "assay", "n_events", "r_all", "n_sig", "r_sig", "n_both_gt0.1", "sign_agree_both"]); sh.to_csv(f"{O}/ceiling_split_half.tsv", sep="\t", index=False)
ca = pd.DataFrame(ca, columns=["cell_line", "target", "shrna_run", "crispr_run", "n_events", "r_all", "n_sig_either", "r_sig_either", "n_sig_both", "sign_agree_both", "same_control"]); ca.to_csv(f"{O}/ceiling_cross_assay.tsv", sep="\t", index=False)
nl = pd.DataFrame(nl, columns=["cell_line", "run_a", "run_b", "same_control", "same_assay", "n_events", "r_all", "n_sig_either", "r_sig_either"]); nl.to_csv(f"{O}/ceiling_null_pairs.tsv", sep="\t", index=False)
print("== split-half within run (rep1 dPSI vs rep2 dPSI) ==\n", sh.groupby(["cell_line", "assay"])[["r_all", "r_sig", "sign_agree_both", "n_sig"]].describe(percentiles=[.1, .5, .9]).round(3).T.to_string())
print("\n== cross-assay same RBP same cell line (shRNA vs CRISPR) ==\n", ca.round(3).to_string(index=False))
print("\n== null: different RBPs, same cell line ==\n", nl.groupby(["cell_line", "same_control"])[["r_all", "r_sig_either"]].describe(percentiles=[.5, .9, .99]).round(3).T.to_string())
