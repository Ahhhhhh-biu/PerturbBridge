"""H4 on the two cohorts added 2026-11-09, with the standards pre-registered for 05_h4_cohort_mapping.py (per-event OLS PSI ~ group +
covariates on rMATS per-sample PSI with coverage >= 10 and >= 70% sample coverage; query = p < 0.05 & |beta| > 0.05, capped at 500;
Pearson r against the dictionary; 200 label permutations -> per-RBP perm_p, p_top = P(max null >= max observed); BH over RBPs).
Cohort 1, GSE182338: FACS-sorted luminal (LEP) and myoepithelial (MEP) human mammary epithelial cells, 51-bp single-end, 120 samples from
48 strains (women); reduction mammoplasty RM (24 strains), contralateral-to-tumour CLTT (12), prophylactic mastectomy PM (8, carriers),
peri-tumour PTT (4). Sorted cells remove the cell-composition confound that dominated the fallopian/breast cohorts. Analysed within
cell type; labels are strain-level so permutations shuffle strains and broadcast to their replicate samples (cluster permutation).
Pre-registered contrasts (each within LEP and within MEP; covariates age z-scored + RNA-seq batch dummies; also unadjusted):
  carrier_vs_RM: pathogenic BRCA1/BRCA2/PALB2 carriers (any tissue type) vs RM; PM_vs_RM: prophylactic mastectomy vs RM;
  adjacent_vs_RM: CLTT + PTT (tumour-bearing breasts, any genotype) vs RM.  Descriptive: LEP_vs_MEP (strong biology; what does it retrieve?).
Cohort 2, GSE289499: adenoma (10) -> carcinoma in situ (8) -> carcinoma (11), Ion Torrent single-end, no normal. Contrasts CIS_vs_adenoma,
carcinoma_vs_adenoma, ordinal trend (0/1/2); adjusted = marker composition scores (epithelial, stroma, immune) + unadjusted.
Outputs: 14_inverse/h4new_<cohort>_<contrast>_<adj>_ranking.tsv, h4new_summary.tsv; stdout -> logs/h4_new_cohorts.txt"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/well/blagden/users/eic789/PerturbBridge/14_inverse"); import retrieval as R
PB = R.PB; O = R.O; rng = np.random.default_rng(0); pd.set_option("display.width", 250); NPERM = 200
P = R.load_dictionary(); info = R.rbp_info(); print("dictionary", P.shape)
def bam_order(d, b): return [p.split("/")[-1].split(".")[0] for p in open(f"{d}/{b}.txt").read().strip().split(",")]
def psi_matrix(d):
    g1, g2 = bam_order(d, "b1"), bam_order(d, "b2"); samples = g1 + g2
    T = R.rmats_psi_table(d); T["sample"] = T["sample"].map({f"g1_{i}": s for i, s in enumerate(g1)} | {f"g2_{i}": s for i, s in enumerate(g2)})
    T.loc[T["cov"] < 10, "psi"] = np.nan; W = T.pivot_table(index="sample", columns="event_id", values="psi", aggfunc="first").reindex(samples)
    return W
def zs(x): x = np.asarray(x, float); return (x - np.nanmean(x)) / (np.nanstd(x) + 1e-12)
def covariates(df, cols):
    out = []
    for v in cols:
        c = df[v]
        if pd.api.types.is_numeric_dtype(c): out.append(zs(c.fillna(c.median()).values)[:, None])
        else: out.append(pd.get_dummies(c.astype(str), drop_first=True).values.astype(float))
    return out
summ = []
def run_contrast(name, W, group, cov_cols, perm_units=None):
    """W: samples x events (NaN allowed), group: 1/0 per sample (float), cov_cols: list of sample x k arrays.
    perm_units: array of cluster ids per sample; labels are permuted across clusters (all samples of a cluster share a label)."""
    W = W.loc[:, W.notna().mean(0) >= 0.7]; Y = W.values; group = np.asarray(group, float)
    def design(grp): return np.hstack([np.ones((len(grp), 1)), grp[:, None]] + cov_cols)
    def query_from(grp):
        b, t, p = R.ols_group_effect(Y, design(grp)); sel = np.where((p < 0.05) & (np.abs(b) > 0.05))[0]
        if len(sel) > 500: sel = sel[np.argsort(p[sel])[:500]]
        return pd.Series(b[sel], index=W.columns[sel]), b, p
    def permute(grp):
        if perm_units is None: return rng.permutation(grp)
        u, idx = np.unique(perm_units, return_inverse=True); lab = np.array([grp[perm_units == k][0] for k in u]); return rng.permutation(lab)[idx]
    q, b_all, p_all = query_from(group); n1 = int(group.sum()); n0 = int(len(group) - group.sum())
    print(f"\n== {name}: samples {len(group)} (group1 {n1}, group0 {n0}; clusters {len(np.unique(perm_units)) if perm_units is not None else len(group)}), events >=70% covered {W.shape[1]}, tested {np.isfinite(p_all).sum()}, p<0.05 {(p_all < 0.05).sum()}, query {len(q)}")
    row = dict(contrast=name, n_samples=len(group), n_group1=n1, n_group0=n0, n_events_tested=int(np.isfinite(p_all).sum()), n_p05=int((p_all < 0.05).sum()), n_query=len(q))
    if len(q) < 50: row["note"] = "too few query events"; summ.append(row); print("   too few query events"); return
    s, n = R.score_against(q, P); rk = R.rank_report(s, info, name)
    null = np.full((NPERM, len(P.columns)), np.nan); nq = []
    for k in range(NPERM):
        qk, _, _ = query_from(permute(group)); nq.append(len(qk))
        if len(qk) >= 50: sk, _ = R.score_against(qk, P); null[k] = sk.reindex(P.columns).values
    obs = s.reindex(P.columns).values; emp_p = np.array([(np.nansum(null[:, j] >= obs[j]) + 1) / (np.isfinite(null[:, j]).sum() + 1) if np.isfinite(obs[j]) else np.nan for j in range(len(obs))])
    top_null = np.nanmax(null, 1); p_top = (np.sum(top_null >= np.nanmax(obs)) + 1) / (np.isfinite(top_null).sum() + 1)
    rk = rk.merge(pd.DataFrame({"rbp": P.columns, "perm_p": emp_p, "null_mean": np.nanmean(null, 0), "null_sd": np.nanstd(null, 0)}), on="rbp")
    rk["z_vs_null"] = (rk.score - rk.null_mean) / rk.null_sd
    m = rk.perm_p.notna(); pv = rk.loc[m, "perm_p"].values; o = np.argsort(pv); bh = np.empty_like(pv); bh[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]; rk.loc[m, "bh_q"] = np.minimum(bh, 1)
    rk.to_csv(f"{O}/h4new_{name}_ranking.tsv", sep="\t", index=False, float_format="%.4f")
    print(f"   top {rk.rbp.iloc[0]} score {rk.score.iloc[0]:.3f} perm_p {rk.perm_p.iloc[0]:.3f} | max observed {np.nanmax(obs):.3f} vs permutation-max null median {np.nanmedian(top_null):.3f} -> p_top {p_top:.3f}; null query size median {np.median(nq):.0f}; perm_p<0.05: {(rk.perm_p < 0.05).sum()} of {s.notna().sum()} (expected {0.05 * s.notna().sum():.0f}); BH q<0.1: {(rk.bh_q < 0.1).sum()}")
    print(rk.head(10)[["rank", "rbp", "score", "perm_p", "bh_q", "z_vs_null", "repro_class"]].round(3).to_string(index=False))
    row.update(n_rbps=int(s.notna().sum()), top_rbp=rk.rbp.iloc[0], top_score=rk.score.iloc[0], top_perm_p=rk.perm_p.iloc[0], p_top_vs_maxnull=p_top, n_rbp_perm_p05=int((rk.perm_p < 0.05).sum()), expected_p05=0.05 * int(s.notna().sum()),
               n_bh_q10=int((rk.bh_q < 0.1).sum()), high_perm_p05=";".join(rk[(rk.repro_class == "high") & (rk.perm_p < 0.05)].rbp)); summ.append(row)

# ---------------- cohort 1: GSE182338 sorted HMEC ----------------
C1 = f"{PB}/05_precancer_cohorts/GSE182338"; d1 = f"{C1}/rmats_carrier_vs_rest"
st = pd.read_csv(f"{C1}/sample_table.tsv", sep="\t").set_index("gsm"); W1 = psi_matrix(d1); st = st.loc[W1.index]
st["carrier"] = st.reported_germline_mutation.isin(["BRCA1", "BRCA2", "PALB2"]); st["age"] = pd.to_numeric(st.age, errors="coerce")
print("\nGSE182338 samples", len(st), "| strains", st.strain_id.nunique(), "| tissue types", st.tissue_type_id.value_counts().to_dict(), "| carriers", int(st.carrier.sum()))
comp1 = pd.read_csv(f"{C1}/cell_composition_scores.tsv", sep="\t").set_index("gsm").reindex(W1.index) if __import__("os").path.exists(f"{C1}/cell_composition_scores.tsv") else None
if comp1 is not None:
    print("marker scores by cell type (mean log2 CPM):"); print(comp1.join(st.cell_type_id).groupby("cell_type_id")[[c for c in ("epithelial", "breast_luminal", "breast_basal", "stroma", "immune", "endothelial") if c in comp1]].mean().round(2).to_string())
for ct in ("LEP", "MEP"):
    sub = st[st.cell_type_id == ct]
    for cname, mask1, mask0 in (("carrier_vs_RM", sub.carrier, sub.tissue_type_id == "RM"), ("PM_vs_RM", sub.tissue_type_id == "PM", sub.tissue_type_id == "RM"),
                                ("adjacent_vs_RM", sub.tissue_type_id.isin(["CLTT", "PTT"]), sub.tissue_type_id == "RM")):
        keep = sub[mask1 | mask0]; grp = mask1[keep.index].astype(float).values
        for adj, cols in (("adjusted", ["age", "batch"]), ("unadjusted", [])):
            cov = covariates(keep, cols); run_contrast(f"GSE182338_{ct}_{cname}_{adj}", W1.loc[keep.index], grp, cov, perm_units=keep.strain_id.values)
        if comp1 is not None:
            cols = [c for c in ("stroma", "immune", "endothelial") if c in comp1]; cov = covariates(keep, ["age", "batch"]) + [zs(comp1.loc[keep.index, c].values)[:, None] for c in cols]
            run_contrast(f"GSE182338_{ct}_{cname}_adjusted_plus_composition", W1.loc[keep.index], grp, cov, perm_units=keep.strain_id.values)
# descriptive: LEP vs MEP within RM (strong, known biology), cluster = strain
rm = st[st.tissue_type_id == "RM"]; run_contrast("GSE182338_RM_LEP_vs_MEP_descriptive", W1.loc[rm.index], (rm.cell_type_id == "LEP").astype(float).values, covariates(rm, ["batch"]), perm_units=rm.strain_id.values)

# ---------------- cohort 2: GSE289499 adenoma -> CIS -> carcinoma ----------------
C2 = f"{PB}/05_precancer_cohorts/GSE289499"; d2 = f"{C2}/rmats_progression"
su = pd.read_csv(f"{C2}/sample_table.tsv", sep="\t").set_index("gsm"); W2 = psi_matrix(d2); su = su.loc[W2.index]
stage = su.tissue.map({"Adenoma": 0, "Carcinoma In Situ": 1, "Carcinoma": 2}); print("\nGSE289499 samples", len(su), su.tissue.value_counts().to_dict())
comp2 = pd.read_csv(f"{C2}/cell_composition_scores.tsv", sep="\t").set_index("gsm").reindex(W2.index) if __import__("os").path.exists(f"{C2}/cell_composition_scores.tsv") else None
if comp2 is not None: print("marker scores by stage:"); print(comp2.join(su.tissue).groupby("tissue")[[c for c in ("epithelial", "stroma", "immune", "endothelial") if c in comp2]].mean().round(2).to_string())
ccols = [c for c in ("epithelial", "stroma", "immune") if comp2 is not None and c in comp2]
for cname, mask1, mask0 in (("CIS_vs_adenoma", stage == 1, stage == 0), ("carcinoma_vs_adenoma", stage == 2, stage == 0), ("advanced_vs_adenoma", stage >= 1, stage == 0)):
    keep = su[mask1 | mask0]; grp = mask1[keep.index].astype(float).values
    run_contrast(f"GSE289499_{cname}_unadjusted", W2.loc[keep.index], grp, [])
    if ccols: run_contrast(f"GSE289499_{cname}_adjusted", W2.loc[keep.index], grp, [zs(comp2.loc[keep.index, c].values)[:, None] for c in ccols])
# ordinal trend: group column = stage (0/1/2), permuted over samples
run_contrast("GSE289499_stage_trend_unadjusted", W2, stage.values.astype(float), [])
if ccols: run_contrast("GSE289499_stage_trend_adjusted", W2, stage.values.astype(float), [zs(comp2[c].values)[:, None] for c in ccols])
S = pd.DataFrame(summ); S.to_csv(f"{O}/h4new_summary.tsv", sep="\t", index=False, float_format="%.4f"); print("\n== summary ==\n", S.round(3).to_string(index=False)); print("DONE")
