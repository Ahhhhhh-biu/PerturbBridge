"""Pre-registered second-stage H4 analyses (STATUS 2026-11-09 16:10). For each cohort:
  A  sample-level inverse mapping: each sample's PSI deviation from the cohort mean is scored against every dictionary RBP (Pearson r over
     >= 50 shared events) -> activity matrix S (samples x RBPs); per RBP OLS  S ~ group + covariates; BH over RBPs; 1,000 group-label permutations
  C  same with the dictionary restricted to each RBP's replicate-consistent ('clean') events (>= 30 shared)
  D  expression layer: STAR ReadsPerGene (unstranded column) -> log2 CPM deviation from cohort mean -> expression dictionary -> same model
  B  (fallopian only) homogeneous subset: premenopausal & batch C; carriers vs non-carriers, BRCA1m vs non, BRCA2m vs non; group-level (as 05) and A
  Positive control for sensitivity: RBP activity scores vs the epithelial composition score (expected to pick mesenchymal/epithelial splicing regulators).
Positive call (pre-registered): BH q < 0.1 AND permutation p < 0.05 AND reproducibility class >= mid.
Outputs: 14_inverse/h4ext_<cohort>_<analysis>.tsv, h4ext_summary.tsv; stdout -> logs/h4_extended.txt"""
import sys, glob, numpy as np, pandas as pd; from scipy import stats; sys.path.insert(0, "/well/blagden/users/eic789/PerturbBridge/14_inverse"); import retrieval as R
PB = R.PB; O = R.O; rng = np.random.default_rng(0); pd.set_option("display.width", 260); NPERM = 1000
P = R.load_dictionary(); info = R.rbp_info()
# clean-event-restricted dictionary
mats = []
for cl in ("HepG2", "K562"):
    g = pd.read_parquet(f"{O}/dict_splicing_rbp_{cl}.parquet"); g = g[g.n_clean >= 1]; mats.append(g.pivot(index="event_id", columns="target", values="dpsi_cc"))
idx = mats[0].index.union(mats[1].index); cols = mats[0].columns.union(mats[1].columns); Pc = pd.DataFrame(np.nanmean(np.stack([m.reindex(index=idx, columns=cols).values for m in mats]), 0), index=idx, columns=cols)
# expression dictionary (base Ensembl ids)
ex = []
for cl in ("HepG2", "K562"):
    e = pd.read_parquet(f"{O}/dict_expression_rbp_{cl}.parquet"); e["gid"] = e.gene_id.str.split(".").str[0]; ex.append(e.pivot_table(index="gid", columns="target", values="dlog2tpm_cc", aggfunc="mean"))
gi = ex[0].index.union(ex[1].index); gc = ex[0].columns.union(ex[1].columns); PE = pd.DataFrame(np.nanmean(np.stack([m.reindex(index=gi, columns=gc).values for m in ex]), 0), index=gi, columns=gc)
print("dictionaries: splicing", P.shape, "clean-restricted", Pc.shape, "expression", PE.shape)
def bh(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); ok = np.isfinite(p); n = ok.sum(); o = np.argsort(p[ok]); r = np.empty(n); r[o] = np.arange(1, n + 1); q[ok] = np.minimum(1, p[ok] * n / r)
    qq = q[ok]; qq = np.minimum.accumulate(qq[o[::-1]])[::-1]; q2 = q.copy(); q2[np.where(ok)[0][o]] = qq; return q2
def activity(W, D, min_shared):
    """W: samples x features (NaN allowed); D: features x RBPs. Returns S samples x RBPs (Pearson r of each sample's deviation from the cohort mean)."""
    dev = W - np.nanmean(W.values, 0); S = np.full((W.shape[0], D.shape[1]), np.nan)
    for i in range(W.shape[0]):
        q = pd.Series(dev.iloc[i].values, index=W.columns).dropna(); s, n = R.score_against(q, D, min_shared=min_shared); S[i] = s.reindex(D.columns).values
    return pd.DataFrame(S, index=W.index, columns=D.columns)
def sample_level(S, group, Xcov, label, cohort):
    X = np.hstack([np.ones((len(group), 1)), group[:, None]] + Xcov); b, t, p = R.ols_group_effect(S.values, X)
    null_t = np.full((NPERM, S.shape[1]), np.nan)
    for k in range(NPERM): null_t[k] = R.ols_group_effect(S.values, np.hstack([np.ones((len(group), 1)), rng.permutation(group)[:, None]] + Xcov))[1]
    perm_p = np.array([(np.nansum(np.abs(null_t[:, j]) >= abs(t[j])) + 1) / (np.isfinite(null_t[:, j]).sum() + 1) if np.isfinite(t[j]) else np.nan for j in range(len(t))])
    d = pd.DataFrame({"rbp": S.columns, "beta": b, "t": t, "p": p, "q_bh": bh(p), "perm_p": perm_p}).set_index("rbp").join(info).reset_index().sort_values("p")
    d["positive_call"] = (d.q_bh < 0.1) & (d.perm_p < 0.05) & d.repro_class.isin(["high", "mid"]); d.to_csv(f"{O}/h4ext_{cohort}_{label}.tsv", sep="\t", index=False, float_format="%.4g")
    maxt_null = np.nanmax(np.abs(null_t), 1); p_fam = (np.sum(maxt_null >= np.nanmax(np.abs(t))) + 1) / (len(maxt_null) + 1)
    print(f"\n-- {cohort} / {label}: n={len(group)} (group1 {int(group.sum())}), RBPs {int(np.isfinite(p).sum())} | min q_bh {np.nanmin(d.q_bh):.3f} | family-wise p of max |t| {p_fam:.3f} | positive calls: {d[d.positive_call].rbp.tolist()}")
    print(d.head(10)[["rbp", "beta", "t", "p", "q_bh", "perm_p", "repro_class"]].round(4).to_string(index=False))
    return dict(cohort=cohort, analysis=label, n=len(group), n_group1=int(group.sum()), n_rbps=int(np.isfinite(p).sum()), top_rbp=d.rbp.iloc[0], top_p=d.p.iloc[0], min_q_bh=float(np.nanmin(d.q_bh)), p_familywise_maxt=p_fam, n_positive_calls=int(d.positive_call.sum()), positive_calls=";".join(d[d.positive_call].rbp))
def group_level(W, group, Xcov, D, label, cohort, nperm=200):
    X = np.hstack([np.ones((len(group), 1)), group[:, None]] + Xcov)
    def query(grp):
        b, t, p = R.ols_group_effect(W.values, np.hstack([np.ones((len(grp), 1)), grp[:, None]] + Xcov)); sel = np.where((p < 0.05) & (np.abs(b) > 0.05))[0]
        if len(sel) > 500: sel = sel[np.argsort(p[sel])[:500]]
        return pd.Series(b[sel], index=W.columns[sel])
    q = query(group); s, n = R.score_against(q, D); obs = s.reindex(D.columns).values; null = np.full((nperm, D.shape[1]), np.nan)
    for k in range(nperm):
        qk = query(rng.permutation(group))
        if len(qk) >= 50: null[k] = R.score_against(qk, D)[0].reindex(D.columns).values
    top_null = np.nanmax(null, 1); p_top = (np.sum(top_null >= np.nanmax(obs)) + 1) / (np.isfinite(top_null).sum() + 1)
    d = pd.DataFrame({"rbp": D.columns, "score": obs, "perm_p": [(np.nansum(null[:, j] >= obs[j]) + 1) / (np.isfinite(null[:, j]).sum() + 1) if np.isfinite(obs[j]) else np.nan for j in range(len(obs))]}).set_index("rbp").join(info).reset_index().sort_values("score", ascending=False)
    d.to_csv(f"{O}/h4ext_{cohort}_{label}.tsv", sep="\t", index=False, float_format="%.4g"); print(f"\n-- {cohort} / {label} (group level): query {len(q)} events, top {d.rbp.iloc[0]} r={d.score.iloc[0]:.3f}, p(top vs permutation max) = {p_top:.3f}")
    return dict(cohort=cohort, analysis=label, n=len(group), n_group1=int(group.sum()), n_rbps=int(np.isfinite(obs).sum()), top_rbp=d.rbp.iloc[0], top_p=p_top, min_q_bh=np.nan, p_familywise_maxt=p_top, n_positive_calls=0, positive_calls="")
def bam_order(d, b): return [p.split("/")[-1].split(".")[0] for p in open(f"{d}/{b}.txt").read().strip().split(",")]
def covmat(df, cols):
    out = []
    for v in cols:
        c = df[v]; out.append(((c - c.mean()) / c.std()).values[:, None] if pd.api.types.is_numeric_dtype(c) else pd.get_dummies(c.astype(str), drop_first=True).values.astype(float))
    return out
def expression_matrix(stardir, samples):
    M = {}
    for s in samples:
        t = pd.read_csv(f"{stardir}/{s}.ReadsPerGene.out.tab", sep="\t", header=None, index_col=0).iloc[4:, 0]; M[s] = t
    M = pd.DataFrame(M).T; M.columns = [c.split(".")[0] for c in M.columns]; cpm = M.div(M.sum(1), axis=0) * 1e6; L = np.log2(cpm + 1); L = L.loc[:, cpm.mean(0) >= 1]; return L
summ = []
# ---------------- fallopian ----------------
C = f"{PB}/05_precancer_cohorts/GSE274628"; d = f"{C}/rmats_carrier_vs_noncarrier"; g1, g2 = bam_order(d, "b1"), bam_order(d, "b2"); samples = g1 + g2
st = pd.read_csv(f"{C}/sample_table.tsv", sep="\t").set_index("gsm").loc[samples]; comp = pd.read_csv(f"{C}/cell_composition_scores.tsv", sep="\t").set_index("gsm").loc[samples]
meta = pd.DataFrame({"carrier": (st.genotype != "Non-BRCAm").astype(float).values, "genotype": st.genotype.values, "batch": st.batch.values, "pre": (st.reproductive_status == "Pre").astype(float).values, "epithelial": comp.epithelial.values}, index=samples)
T = R.rmats_psi_table(d); T["sample"] = T["sample"].map({f"g1_{i}": s for i, s in enumerate(g1)} | {f"g2_{i}": s for i, s in enumerate(g2)}); T.loc[T["cov"] < 10, "psi"] = np.nan
W = T.pivot_table(index="sample", columns="event_id", values="psi", aggfunc="first").reindex(samples); W = W.loc[:, W.notna().mean(0) >= 0.7]; print("\n== GSE274628 fallopian:", W.shape)
S = activity(W, P, 50); Sc = activity(W, Pc, 30); L = expression_matrix(f"{C}/star", samples); SE = activity(L, PE, 200); print("activity matrices:", S.shape, Sc.shape, "expression genes", L.shape[1], SE.shape)
# positive control: activity vs epithelial score
r_epi = pd.Series({c: stats.pearsonr(S[c].values[np.isfinite(S[c].values)], meta.epithelial.values[np.isfinite(S[c].values)])[0] for c in S.columns}).sort_values(); print("\nsensitivity control, splicing activity vs epithelial score: most negative", r_epi.head(6).round(2).to_dict(), "| most positive", r_epi.tail(6).round(2).to_dict())
cov = covmat(meta, ["batch", "pre", "epithelial"])
summ.append(sample_level(S, meta.carrier.values, cov, "A_sample_level", "fallopian")); summ.append(sample_level(Sc, meta.carrier.values, cov, "C_clean_dictionary", "fallopian")); summ.append(sample_level(SE, meta.carrier.values, cov, "D_expression", "fallopian"))
sub = meta[(meta.pre == 1) & (meta.batch == "C")]; print("\nhomogeneous subset (premenopausal, batch C):", sub.genotype.value_counts().to_dict())
cov_s = covmat(sub, ["epithelial"]); summ.append(group_level(W.loc[sub.index], sub.carrier.values, cov_s, P, "B_subset_group_level", "fallopian")); summ.append(sample_level(S.loc[sub.index], sub.carrier.values, cov_s, "B_subset_sample_level", "fallopian"))
for gt in ("BRCA1m", "BRCA2m"):
    m = sub[sub.genotype.isin([gt, "Non-BRCAm"])]; summ.append(sample_level(S.loc[m.index], (m.genotype == gt).astype(float).values, covmat(m, ["epithelial"]), f"B_subset_{gt}_sample_level", "fallopian"))
# ---------------- breast ----------------
C = f"{PB}/05_precancer_cohorts/GSE166044"; d = f"{C}/rmats_susceptible_vs_healthy"; g1, g2 = bam_order(d, "b1"), bam_order(d, "b2"); samples = g1 + g2
comp = pd.read_csv(f"{C}/cell_composition_scores.tsv", sep="\t").set_index("gsm").loc[samples]; grp = np.array([1.0] * len(g1) + [0.0] * len(g2))
T = R.rmats_psi_table(d); T["sample"] = T["sample"].map({f"g1_{i}": s for i, s in enumerate(g1)} | {f"g2_{i}": s for i, s in enumerate(g2)}); T.loc[T["cov"] < 10, "psi"] = np.nan
W = T.pivot_table(index="sample", columns="event_id", values="psi", aggfunc="first").reindex(samples); W = W.loc[:, W.notna().mean(0) >= 0.7]; print("\n== GSE166044 breast:", W.shape)
S = activity(W, P, 50); Sc = activity(W, Pc, 30); L = expression_matrix(f"{C}/star", samples); SE = activity(L, PE, 200); cov = covmat(comp, ["epithelial"])
summ.append(sample_level(S, grp, cov, "A_sample_level", "breast")); summ.append(sample_level(Sc, grp, cov, "C_clean_dictionary", "breast")); summ.append(sample_level(SE, grp, cov, "D_expression", "breast"))
Sm = pd.DataFrame(summ); Sm.to_csv(f"{O}/h4ext_summary.tsv", sep="\t", index=False, float_format="%.4g"); print("\n== summary ==\n", Sm.round(4).to_string(index=False)); print("DONE")
