"""Block 3: H4 inverse mapping of the two precancer cohorts onto the RBP dictionary.
Cohorts (rMATS outputs already on disk): GSE274628 fallopian tube, BRCA carrier (41) vs non-carrier (53), covariates batch + epithelial score;
GSE166044 normal breast, later-diagnosed 'Susceptible' (15) vs healthy (15), covariate epithelial score.
Per event: OLS PSI ~ group + covariates on samples with coverage >= 10, events covered in >= 70% of samples. Query = events with
p < 0.05 and |beta| > 0.05 (capped at the 500 smallest p); score every dictionary RBP by Pearson r with its mean dpsi_cc over the
shared events (>= 50). Null 1: 200 permutations of the group label (recompute betas, re-select, re-score) -> empirical p per RBP and
for the top score. Null 2: random RBP (uniform ranks). Also an unadjusted version (no covariates) to show the composition effect.
Outputs: 14_inverse/h4_<cohort>_ranking.tsv, h4_<cohort>_query_events.tsv, h4_summary.tsv; stdout -> logs/h4_cohorts.txt"""
import sys, numpy as np, pandas as pd; sys.path.insert(0, "/well/blagden/users/eic789/PerturbBridge/14_inverse"); import retrieval as R
PB = R.PB; O = R.O; rng = np.random.default_rng(0); pd.set_option("display.width", 250); NPERM = 200
P = R.load_dictionary(); info = R.rbp_info(); print("dictionary", P.shape)
def bam_order(d, b):
    return [p.split("/")[-1].split(".")[0] for p in open(f"{d}/{b}.txt").read().strip().split(",")]
COH = {"GSE274628_fallopian": dict(dir=f"{PB}/05_precancer_cohorts/GSE274628/rmats_carrier_vs_noncarrier", comp=f"{PB}/05_precancer_cohorts/GSE274628/cell_composition_scores.tsv", covars=["batch", "epithelial"]),
       "GSE166044_breast": dict(dir=f"{PB}/05_precancer_cohorts/GSE166044/rmats_susceptible_vs_healthy", comp=f"{PB}/05_precancer_cohorts/GSE166044/cell_composition_scores.tsv", covars=["epithelial"])}
summ = []
for name, c in COH.items():
    g1, g2 = bam_order(c["dir"], "b1"), bam_order(c["dir"], "b2"); samples = g1 + g2; group = np.array([1] * len(g1) + [0] * len(g2), float)
    comp = pd.read_csv(c["comp"], sep="\t").set_index("gsm").loc[samples]
    T = R.rmats_psi_table(c["dir"]); T["sample"] = T["sample"].map({f"g1_{i}": s for i, s in enumerate(g1)} | {f"g2_{i}": s for i, s in enumerate(g2)})
    T.loc[T["cov"] < 10, "psi"] = np.nan; W = T.pivot_table(index="sample", columns="event_id", values="psi", aggfunc="first").reindex(samples)
    W = W.loc[:, W.notna().mean(0) >= 0.7]; print(f"\n== {name}: samples {len(samples)} (group1 {len(g1)}), events covered >=70%: {W.shape[1]}")
    Y = W.values
    for adj in ("adjusted", "unadjusted"):
        cov_cols = []
        if adj == "adjusted":
            for v in c["covars"]:
                col = comp[v]; cov_cols.append(((col - col.mean()) / col.std()).values[:, None] if pd.api.types.is_numeric_dtype(col) else pd.get_dummies(col.astype(str), drop_first=True).values.astype(float))
        def design(grp): return np.hstack([np.ones((len(grp), 1)), grp[:, None]] + cov_cols)
        def query_from(grp):
            b, t, p = R.ols_group_effect(Y, design(grp)); sel = np.where((p < 0.05) & (np.abs(b) > 0.05))[0]
            if len(sel) > 500: sel = sel[np.argsort(p[sel])[:500]]
            return pd.Series(b[sel], index=W.columns[sel]), b, p
        q, b_all, p_all = query_from(group); s, n = R.score_against(q, P); rk = R.rank_report(s, info, f"{name}_{adj}")
        print(f"-- {adj}: query events {len(q)} (of {np.isfinite(p_all).sum()} tested; p<0.05: {(p_all < 0.05).sum()}), RBPs scored {s.notna().sum()}")
        # null 1: group-label permutations
        null = np.full((NPERM, len(P.columns)), np.nan); nq = []
        for k in range(NPERM):
            qk, _, _ = query_from(rng.permutation(group)); nq.append(len(qk))
            if len(qk) >= 50: sk, _ = R.score_against(qk, P); null[k] = sk.reindex(P.columns).values
        obs = s.reindex(P.columns).values; emp_p = np.array([(np.nansum(null[:, j] >= obs[j]) + 1) / (np.isfinite(null[:, j]).sum() + 1) if np.isfinite(obs[j]) else np.nan for j in range(len(obs))])
        top_null = np.nanmax(null, 1); p_top = (np.sum(top_null >= np.nanmax(obs)) + 1) / (np.isfinite(top_null).sum() + 1)
        rk = rk.merge(pd.DataFrame({"rbp": P.columns, "perm_p": emp_p, "null_mean": np.nanmean(null, 0), "null_sd": np.nanstd(null, 0)}), on="rbp")
        rk["z_vs_null"] = (rk.score - rk.null_mean) / rk.null_sd; rk.to_csv(f"{O}/h4_{name}_{adj}_ranking.tsv", sep="\t", index=False, float_format="%.4f")
        if adj == "adjusted": pd.DataFrame({"event_id": W.columns, "beta": b_all, "p": p_all}).dropna().to_csv(f"{O}/h4_{name}_query_events.tsv", sep="\t", index=False, float_format="%.5f")
        print(f"   top score {np.nanmax(obs):.3f} vs permutation-max null median {np.nanmedian(top_null):.3f} -> p_top = {p_top:.3f}; median query size under null {np.median(nq):.0f}")
        print(rk.head(15).round(3).to_string(index=False))
        print("   high-reproducibility RBPs with perm_p < 0.05:", rk[(rk.repro_class == "high") & (rk.perm_p < 0.05)].rbp.tolist())
        summ.append(dict(cohort=name, adjustment=adj, n_samples=len(samples), n_events_tested=int(np.isfinite(p_all).sum()), n_query=len(q), n_rbps=int(s.notna().sum()), top_rbp=rk.rbp.iloc[0], top_score=rk.score.iloc[0], top_perm_p=rk.perm_p.iloc[0], p_top_vs_maxnull=p_top,
                         n_rbp_perm_p05=int((rk.perm_p < 0.05).sum()), n_high_perm_p05=int(((rk.repro_class == "high") & (rk.perm_p < 0.05)).sum()), expected_p05=0.05 * int(s.notna().sum())))
pd.DataFrame(summ).to_csv(f"{O}/h4_summary.tsv", sep="\t", index=False, float_format="%.4f"); print("\n", pd.DataFrame(summ).round(3).to_string(index=False)); print("DONE")
