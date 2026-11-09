"""Robustness checks for the second-stage H4 positives (not pre-registered as calls; reported in full):
 R1 full adjustment: all ten cell-composition scores + age (+ batch, menopause for fallopian) instead of epithelial only
 R2 one axis or many: PCA of the activity matrix restricted to the positive-call RBPs; variance explained by PC1; correlation of PC1 with group
 R3 composition dependence: correlation of each positive RBP's activity with each composition score
 R4 direction: sign of beta (positive = carriers resemble knockdown of the RBP; negative = carriers resemble the opposite)
Also saves the activity matrices (h4ext_activity_<cohort>_<layer>.parquet) for reuse.
Outputs: 14_inverse/h4rob_<cohort>_<layer>_fulladj.tsv, h4rob_summary.tsv; stdout -> logs/h4_robustness.txt"""
import sys, numpy as np, pandas as pd; from scipy import stats; sys.path.insert(0, "/well/blagden/users/eic789/PerturbBridge/14_inverse"); import retrieval as R
PB = R.PB; O = R.O; rng = np.random.default_rng(1); pd.set_option("display.width", 260); NPERM = 1000
P = R.load_dictionary(); info = R.rbp_info()
ex = []
for cl in ("HepG2", "K562"):
    e = pd.read_parquet(f"{O}/dict_expression_rbp_{cl}.parquet"); e["gid"] = e.gene_id.str.split(".").str[0]; ex.append(e.pivot_table(index="gid", columns="target", values="dlog2tpm_cc", aggfunc="mean"))
gi = ex[0].index.union(ex[1].index); gc = ex[0].columns.union(ex[1].columns); PE = pd.DataFrame(np.nanmean(np.stack([m.reindex(index=gi, columns=gc).values for m in ex]), 0), index=gi, columns=gc)
def bh(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); ok = np.isfinite(p); n = ok.sum(); o = np.argsort(p[ok]); r = np.empty(n); r[o] = np.arange(1, n + 1); qq = np.minimum(1, p[ok] * n / r); qq = np.minimum.accumulate(qq[o[::-1]])[::-1]; q[np.where(ok)[0][o]] = qq; return q
def activity(W, D, min_shared):
    dev = W - np.nanmean(W.values, 0); S = np.full((W.shape[0], D.shape[1]), np.nan)
    for i in range(W.shape[0]):
        q = pd.Series(dev.iloc[i].values, index=W.columns).dropna(); s, n = R.score_against(q, D, min_shared=min_shared); S[i] = s.reindex(D.columns).values
    return pd.DataFrame(S, index=W.index, columns=D.columns)
def bam_order(d, b): return [p.split("/")[-1].split(".")[0] for p in open(f"{d}/{b}.txt").read().strip().split(",")]
def covmat(df, cols):
    out = []
    for v in cols:
        c = df[v]; out.append(((c - c.mean()) / (c.std() + 1e-9)).values[:, None] if pd.api.types.is_numeric_dtype(c) else pd.get_dummies(c.astype(str), drop_first=True).values.astype(float))
    return out
def expression_matrix(stardir, samples):
    M = pd.DataFrame({s: pd.read_csv(f"{stardir}/{s}.ReadsPerGene.out.tab", sep="\t", header=None, index_col=0).iloc[4:, 0] for s in samples}).T; M.columns = [c.split(".")[0] for c in M.columns]
    cpm = M.div(M.sum(1), axis=0) * 1e6; L = np.log2(cpm + 1); return L.loc[:, cpm.mean(0) >= 1]
def test(S, group, Xcov, label, cohort, prev_calls):
    X = np.hstack([np.ones((len(group), 1)), group[:, None]] + Xcov); b, t, p = R.ols_group_effect(S.values, X)
    null = np.full((NPERM, S.shape[1]), np.nan)
    for k in range(NPERM): null[k] = R.ols_group_effect(S.values, np.hstack([np.ones((len(group), 1)), rng.permutation(group)[:, None]] + Xcov))[1]
    perm_p = np.array([(np.nansum(np.abs(null[:, j]) >= abs(t[j])) + 1) / (np.isfinite(null[:, j]).sum() + 1) if np.isfinite(t[j]) else np.nan for j in range(len(t))])
    d = pd.DataFrame({"rbp": S.columns, "beta": b, "t": t, "p": p, "q_bh": bh(p), "perm_p": perm_p}).set_index("rbp").join(info).reset_index(); d["positive_call"] = (d.q_bh < 0.1) & (d.perm_p < 0.05) & d.repro_class.isin(["high", "mid"]); d["was_positive_before"] = d.rbp.isin(prev_calls)
    d = d.sort_values("p"); d.to_csv(f"{O}/h4rob_{cohort}_{label}_fulladj.tsv", sep="\t", index=False, float_format="%.4g"); pf = (np.sum(np.nanmax(np.abs(null), 1) >= np.nanmax(np.abs(t))) + 1) / (NPERM + 1)
    surv = d[d.was_positive_before & d.positive_call].rbp.tolist(); lost = d[d.was_positive_before & ~d.positive_call].rbp.tolist()
    print(f"\n-- {cohort} / {label} FULL adjustment (p = {X.shape[1]} columns): family-wise p {pf:.3f} | positive calls {int(d.positive_call.sum())} | previous calls surviving {len(surv)}/{len(prev_calls)}: {surv} | lost: {lost[:15]}")
    print(d.head(8)[["rbp", "beta", "t", "q_bh", "perm_p", "repro_class"]].round(4).to_string(index=False))
    return dict(cohort=cohort, layer=label, n_cov=X.shape[1] - 2, p_familywise=pf, n_calls_full=int(d.positive_call.sum()), n_prev=len(prev_calls), n_surviving=len(surv), surviving=";".join(surv))
def axis(S, group, calls, comp, label):
    A = S[[c for c in calls if c in S.columns]].dropna(axis=1); A = A.loc[:, A.notna().all()] if A.notna().all().any() else A
    Z = (A - A.mean()) / A.std(); U, s, Vt = np.linalg.svd(Z.fillna(0).values, full_matrices=False); ve = s ** 2 / (s ** 2).sum(); pc1 = U[:, 0] * s[0]
    r_g = stats.pearsonr(pc1, group)[0]; rc = {c: round(stats.pearsonr(pc1, comp[c].values)[0], 2) for c in comp.columns if pd.api.types.is_numeric_dtype(comp[c])}
    print(f"   [{label}] {A.shape[1]} positive RBPs: PC1 explains {ve[0]:.0%} of their activity variance (PC2 {ve[1]:.0%}); r(PC1, group) = {r_g:+.2f}; r(PC1, composition scores) = {rc}")
    cm = {c: round(np.nanmean([abs(stats.pearsonr(S[r].values, comp[c].values)[0]) for r in A.columns]), 2) for c in comp.columns if pd.api.types.is_numeric_dtype(comp[c])}; print(f"   mean |r| of positive RBPs' activity with each composition score: {cm}")
summ = []
# fallopian
C = f"{PB}/05_precancer_cohorts/GSE274628"; d = f"{C}/rmats_carrier_vs_noncarrier"; g1, g2 = bam_order(d, "b1"), bam_order(d, "b2"); samples = g1 + g2
st = pd.read_csv(f"{C}/sample_table.tsv", sep="\t").set_index("gsm").loc[samples]; comp = pd.read_csv(f"{C}/cell_composition_scores.tsv", sep="\t").set_index("gsm").loc[samples]; compcols = ["epithelial", "ft_secretory", "ft_ciliated", "breast_luminal", "breast_basal", "stroma", "smooth_muscle", "immune", "endothelial", "adipose"]
meta = pd.concat([comp[compcols], pd.DataFrame({"carrier": (st.genotype != "Non-BRCAm").astype(float).values, "batch": st.batch.values, "pre": (st.reproductive_status == "Pre").astype(float).values, "age": pd.to_numeric(st.age_at_time_of_surgery, errors="coerce").fillna(pd.to_numeric(st.age_at_time_of_surgery, errors="coerce").median()).values}, index=samples)], axis=1)
T = R.rmats_psi_table(d); T["sample"] = T["sample"].map({f"g1_{i}": s for i, s in enumerate(g1)} | {f"g2_{i}": s for i, s in enumerate(g2)}); T.loc[T["cov"] < 10, "psi"] = np.nan
W = T.pivot_table(index="sample", columns="event_id", values="psi", aggfunc="first").reindex(samples); W = W.loc[:, W.notna().mean(0) >= 0.7]
S = activity(W, P, 50); SE = activity(expression_matrix(f"{C}/star", samples), PE, 200); S.to_parquet(f"{O}/h4ext_activity_fallopian_splicing.parquet"); SE.to_parquet(f"{O}/h4ext_activity_fallopian_expression.parquet")
full = covmat(meta, ["batch", "pre", "age"] + compcols); print("== fallopian: full covariates =", 3 + len(compcols), "columns ->", sum(x.shape[1] for x in full))
for layer, M in (("A_sample_level", S), ("D_expression", SE)):
    prev = pd.read_csv(f"{O}/h4ext_fallopian_{layer}.tsv", sep="\t"); calls = prev[prev.positive_call].rbp.tolist()
    summ.append(test(M, meta.carrier.values, full, layer, "fallopian", calls)); axis(M, meta.carrier.values, calls, comp[compcols], f"fallopian {layer}")
# breast
C = f"{PB}/05_precancer_cohorts/GSE166044"; d = f"{C}/rmats_susceptible_vs_healthy"; g1, g2 = bam_order(d, "b1"), bam_order(d, "b2"); samples = g1 + g2
comp = pd.read_csv(f"{C}/cell_composition_scores.tsv", sep="\t").set_index("gsm").loc[samples]; st = pd.read_csv(f"{C}/sample_table.tsv", sep="\t").set_index("gsm").loc[samples]; grp = np.array([1.0] * len(g1) + [0.0] * len(g2))
meta = comp[compcols].copy(); meta["age"] = pd.to_numeric(st.age, errors="coerce").fillna(50).values; meta["bmi"] = pd.to_numeric(st.body_mass_index, errors="coerce").fillna(25).values
T = R.rmats_psi_table(d); T["sample"] = T["sample"].map({f"g1_{i}": s for i, s in enumerate(g1)} | {f"g2_{i}": s for i, s in enumerate(g2)}); T.loc[T["cov"] < 10, "psi"] = np.nan
W = T.pivot_table(index="sample", columns="event_id", values="psi", aggfunc="first").reindex(samples); W = W.loc[:, W.notna().mean(0) >= 0.7]
S = activity(W, P, 50); SE = activity(expression_matrix(f"{C}/star", samples), PE, 200); S.to_parquet(f"{O}/h4ext_activity_breast_splicing.parquet"); SE.to_parquet(f"{O}/h4ext_activity_breast_expression.parquet")
full = covmat(meta, ["age", "bmi"] + compcols[:1] + ["breast_luminal", "breast_basal", "stroma", "immune", "endothelial", "adipose"]); print("\n== breast: full covariates ->", sum(x.shape[1] for x in full), "columns (n = 30, so kept to 9)")
for layer, M in (("A_sample_level", S), ("D_expression", SE)):
    prev = pd.read_csv(f"{O}/h4ext_breast_{layer}.tsv", sep="\t"); calls = prev[prev.positive_call].rbp.tolist()
    summ.append(test(M, grp, full, layer, "breast", calls))
    if calls: axis(M, grp, calls, comp[compcols], f"breast {layer}")
pd.DataFrame(summ).to_csv(f"{O}/h4rob_summary.tsv", sep="\t", index=False, float_format="%.4g"); print("\n", pd.DataFrame(summ).to_string(index=False)); print("DONE")
