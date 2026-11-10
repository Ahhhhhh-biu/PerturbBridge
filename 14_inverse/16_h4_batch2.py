"""Inverse mapping on the second batch of cohorts (plan and success criteria: 14_inverse/H4_BATCH2_PLAN.md, written before the results).
Usage: python 16_h4_batch2.py <GSE>   (GSE220518 | GSE193946 | GSE79209 | GSE227919 | GSE202048); reads <cohort>/rmats_all/ (rMATS post over all samples)
Outputs: 14_inverse/h4b2_<contrast>_ranking.tsv, h4b2_summary_<GSE>.tsv, 05_precancer_cohorts/<GSE>/marker_scores.tsv; stdout -> logs/h4_batch2_<GSE>.txt"""
import sys, re, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "/well/blagden/users/eic789/PerturbBridge/14_inverse"); import retrieval as R
PB = R.PB; O = R.O; rng = np.random.default_rng(0); pd.set_option("display.width", 250); NPERM = 200; gse = sys.argv[1]; d0 = f"{PB}/05_precancer_cohorts/{gse}"; d = f"{d0}/rmats_all"
P = R.load_dictionary(); info = R.rbp_info(); print("dictionary", P.shape, flush=True); GTF = "/well/blagden/users/eic789/LncRNASPLICE/gencode.v47.primary_assembly.annotation.gtf"
SF3 = ["SF3B1", "SF3B3", "SF3B4", "SF3B2", "SF3A1", "SF3A2", "SF3A3"]
def bam_order(b): return [p.split("/")[-1].split(".")[0] for p in open(f"{d}/{b}.txt").read().strip().split(",")]
def psi_matrix():
    g1, g2 = bam_order("b1"), bam_order("b2"); samples = g1 + g2; T = R.rmats_psi_table(d); T["sample"] = T["sample"].map({f"g1_{i}": s for i, s in enumerate(g1)} | {f"g2_{i}": s for i, s in enumerate(g2)})
    T.loc[T["cov"] < 10, "psi"] = np.nan; return T.pivot_table(index="sample", columns="event_id", values="psi", aggfunc="first").reindex(samples)
def marker_scores(samples, MARK):
    sym = {}
    for l in open(GTF):
        if l.startswith("#"): continue
        f = l.split("\t")
        if f[2] == "gene": sym[re.search('gene_id "([^"]+)"', f[8]).group(1)] = re.search('gene_name "([^"]+)"', f[8]).group(1)
    nf = pd.read_csv(f"{d0}/star/{samples[0]}.ReadsPerGene.out.tab", sep="\t", header=None, nrows=4, index_col=0).loc["N_noFeature"]; col = 2 if nf[2] < nf[3] / 3 else 3 if nf[3] < nf[2] / 3 else 1
    C = pd.DataFrame({g: pd.read_csv(f"{d0}/star/{g}.ReadsPerGene.out.tab", sep="\t", header=None, index_col=0, skiprows=4)[col] for g in samples}); C.index = C.index.map(sym); C = C.groupby(level=0).sum(); L = np.log2(C / C.sum() * 1e6 + 1)
    S = pd.DataFrame({k: L.reindex(v).dropna().mean() for k, v in MARK.items()}); S.to_csv(f"{d0}/marker_scores.tsv", sep="\t", float_format="%.3f"); print("count column", col, "| marker scores", S.shape, flush=True); return S
def zs(x): x = np.asarray(x, float); return (x - np.nanmean(x)) / (np.nanstd(x) + 1e-12)
def covariates(df, cols):
    out = []
    for v in cols:
        c = df[v]
        if pd.api.types.is_numeric_dtype(c): out.append(zs(c.fillna(c.median()).values)[:, None])
        else: out.append(pd.get_dummies(c.astype(str), drop_first=True).values.astype(float))
    return out
summ = []
def finish(name, q, query_null, row, focus):
    """Score the observed query, build the permutation null with query_null(), write the ranking and the summary row."""
    if len(q) < 50: row["note"] = "too few query events"; summ.append(row); print("   too few query events", flush=True); return
    s, n = R.score_against(q, P); rk = R.rank_report(s, info, name); null = np.full((NPERM, len(P.columns)), np.nan); nq = []
    for k in range(NPERM):
        qk = query_null(); nq.append(len(qk))
        if len(qk) >= 50: sk, _ = R.score_against(qk, P); null[k] = sk.reindex(P.columns).values
    obs = s.reindex(P.columns).values; emp_p = np.array([(np.nansum(null[:, j] >= obs[j]) + 1) / (np.isfinite(null[:, j]).sum() + 1) if np.isfinite(obs[j]) else np.nan for j in range(len(obs))])
    top_null = np.nanmax(null, 1); p_top = (np.sum(top_null >= np.nanmax(obs)) + 1) / (np.isfinite(top_null).sum() + 1)
    rk = rk.merge(pd.DataFrame({"rbp": P.columns, "perm_p": emp_p, "null_mean": np.nanmean(null, 0), "null_sd": np.nanstd(null, 0)}), on="rbp"); rk["z_vs_null"] = (rk.score - rk.null_mean) / rk.null_sd
    m = rk.perm_p.notna(); pv = rk.loc[m, "perm_p"].values; o = np.argsort(pv); bh = np.empty_like(pv); bh[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]; rk.loc[m, "bh_q"] = np.minimum(bh, 1)
    rk.to_csv(f"{O}/h4b2_{name}_ranking.tsv", sep="\t", index=False, float_format="%.4f"); nr = int(s.notna().sum()); top10 = list(rk.rbp.head(10)); hi = rk[(rk.repro_class == "high") & (rk.perm_p < 0.05)].rbp.tolist()
    print(f"   top {rk.rbp.iloc[0]} score {rk.score.iloc[0]:.3f} perm_p {rk.perm_p.iloc[0]:.3f} | max observed {np.nanmax(obs):.3f} vs permutation-max null median {np.nanmedian(top_null):.3f} (95th pct {np.nanpercentile(top_null, 95):.3f}) -> p_top {p_top:.3f}; null query size median {np.median(nq):.0f}; perm_p<0.05: {(rk.perm_p < 0.05).sum()} of {nr} (expected {0.05 * nr:.0f}); BH q<0.1: {(rk.bh_q < 0.1).sum()}")
    print(rk.head(10)[["rank", "rbp", "score", "perm_p", "bh_q", "z_vs_null", "repro_class"]].round(3).to_string(index=False))
    fr = {g: (int(rk.loc[rk.rbp == g, "rank"].iloc[0]), float(rk.loc[rk.rbp == g, "perm_p"].iloc[0])) for g in focus if (rk.rbp == g).any()}
    if fr: print("   focus: " + "; ".join(f"{g} rank {r_} of {nr} (perm p {p_:.3f})" for g, (r_, p_) in fr.items()) + f" | SF3A/SF3B members in the top ten: {sum(g in SF3 for g in top10)}", flush=True)
    row.update(n_rbps=nr, top_rbp=rk.rbp.iloc[0], top_score=rk.score.iloc[0], top_perm_p=rk.perm_p.iloc[0], p_top_vs_maxnull=p_top, null_top_median=np.nanmedian(top_null), null_top_p95=np.nanpercentile(top_null, 95), n_rbp_perm_p05=int((rk.perm_p < 0.05).sum()), expected_p05=0.05 * nr,
               n_bh_q10=int((rk.bh_q < 0.1).sum()), high_perm_p05=";".join(hi), focus=";".join(f"{g}:{r_}" for g, (r_, p_) in fr.items()), n_sf3_top10=sum(g in SF3 for g in top10)); summ.append(row)
def select(b, p):
    sel = np.where((p < 0.05) & (np.abs(b) > 0.05))[0]
    return sel[np.argsort(p[sel])[:500]] if len(sel) > 500 else sel
def run_contrast(name, W, group, cov_cols, focus=()):
    """W: samples x events (NaN allowed); group: 1/0 or ordinal per sample; cov_cols: list of sample x k arrays."""
    W = W.loc[:, W.notna().mean(axis=0) >= 0.7]; Y = W.values; group = np.asarray(group, float)
    def query_from(grp): b, t, p = R.ols_group_effect(Y, np.hstack([np.ones((len(grp), 1)), grp[:, None]] + cov_cols)); sel = select(b, p); return pd.Series(b[sel], index=W.columns[sel]), p
    q, p_all = query_from(group); vals = np.unique(group); n1 = int((group == vals.max()).sum()); n0 = int((group == vals.min()).sum())
    print(f"\n== {name}: samples {len(group)} (highest level {n1}, lowest {n0}; levels {len(vals)}), events >=70% covered {W.shape[1]}, tested {np.isfinite(p_all).sum()}, p<0.05 {(p_all < 0.05).sum()}, query {len(q)}", flush=True)
    finish(name, q, lambda: query_from(rng.permutation(group))[0], dict(contrast=name, n_samples=len(group), n_group1=n1, n_group0=n0, n_events_tested=int(np.isfinite(p_all).sum()), n_p05=int((p_all < 0.05).sum()), n_query=len(q)), focus)
def run_paired(name, W, case, ctrl, dcov=None, focus=()):
    """Paired design: D = PSI(case) - PSI(control) per pair; effect = mean of D (optionally adjusted for paired differences in covariates, no intercept);
    null = random sign flips of pairs."""
    Dm = W.loc[case].values - W.loc[ctrl].values; keep = np.isfinite(Dm).mean(0) >= 0.7; Dm = Dm[:, keep]; cols = W.columns[keep]; npair = len(case)
    def query_from(s):
        if dcov is None:
            Ds = Dm * s[:, None]; n = np.isfinite(Ds).sum(0); mu = np.nanmean(Ds, 0); sd = np.nanstd(Ds, 0, ddof=1); t = mu / (sd / np.sqrt(n) + 1e-12); p = np.where(sd > 1e-10, 2 * stats.t.sf(np.abs(t), n - 1), np.nan); b = mu
        else: b, t, p = R.ols_group_effect(Dm, np.hstack([dcov[:, :1], s[:, None], dcov[:, 1:]]))
        sel = select(b, np.where(np.isfinite(p), p, 1)); return pd.Series(b[sel], index=cols[sel]), p
    q, p_all = query_from(np.ones(npair)); print(f"\n== {name}: pairs {npair}, events >=70% covered {keep.sum()}, tested {np.isfinite(p_all).sum()}, p<0.05 {(p_all < 0.05).sum()}, query {len(q)}", flush=True)
    finish(name, q, lambda: query_from(rng.choice([-1.0, 1.0], npair))[0], dict(contrast=name, n_samples=2 * npair, n_group1=npair, n_group0=npair, n_events_tested=int(np.isfinite(p_all).sum()), n_p05=int((p_all < 0.05).sum()), n_query=len(q)), focus)
st = pd.read_csv(f"{d0}/sample_table.tsv", sep="\t").set_index("gsm"); W = psi_matrix(); st = st.reindex(W.index); print(gse, "samples", len(st), "| PSI matrix", W.shape, flush=True)
EPI = {"epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"], "squamous": ["KRT5", "KRT14", "KRT13", "KRT4"], "stroma": ["COL1A1", "COL1A2", "DCN", "LUM"], "immune": ["PTPRC"], "endothelial": ["PECAM1", "VWF", "CDH5"], "smooth_muscle": ["ACTA2", "MYH11", "DES", "CNN1"], "ciliated": ["FOXJ1", "CAPS", "TPPP3", "PIFO"]}
def sub(mask): return W[mask.values], st[mask.values]
if gse == "GSE220518":
    S = marker_scores(list(st.index), {"stem": ["CD34", "AVP", "HLF", "CRHBP"], "erythroid": ["HBB", "GYPA", "KLF1", "GATA1"], "myeloid": ["MPO", "ELANE", "AZU1", "PRTN3"], "lymphoid": ["CD79A", "VPREB1", "DNTT", "IGLL1"], "megakaryocytic": ["PF4", "ITGA2B", "GP9"]}); st = st.join(S); lin = list(S.columns)
    g = st.genotype.fillna(""); st["sf3b1"] = g.str.contains("SF3B1 MUT").astype(float); st["srsf2"] = g.str.contains("SRSF2 MUT").astype(float); st["u2af1"] = g.str.contains("U2AF1 MUT").astype(float); st["nmut"] = st.sf3b1 + st.srsf2 + st.u2af1; mds = st.disease == "LR-MDS"
    print("MDS", int(mds.sum()), "| SF3B1 only", int((mds & (st.sf3b1 == 1) & (st.nmut == 1)).sum()), "SRSF2 only", int((mds & (st.srsf2 == 1) & (st.nmut == 1)).sum()), "U2AF1 only", int((mds & (st.u2af1 == 1) & (st.nmut == 1)).sum()), "no SF mutation", int((mds & (st.nmut == 0)).sum()), "| progression", dict(st.progression.value_counts()))
    print("lineage scores by group:\n", st[mds].assign(grp=np.where(st[mds].nmut == 0, "WT", np.where(st[mds].sf3b1 == 1, "SF3B1", np.where(st[mds].srsf2 == 1, "SRSF2", "U2AF1")))).groupby("grp")[lin].mean().round(2).to_string())
    for key, col, foc in (("A1_SF3B1mut_vs_WT", "sf3b1", ("SF3B1", "SF3B3", "SF3A1", "U2AF1", "U2AF2")), ("A2_U2AF1mut_vs_WT", "u2af1", ("U2AF1", "U2AF2", "SF3B1")), ("A3_SRSF2mut_vs_WT", "srsf2", ("SRSF1", "SRSF3", "SRSF5", "SRSF7", "SF3B1"))):
        m = mds & (((st[col] == 1) & (st.nmut == 1)) | (st.nmut == 0)); Wm, sm = sub(m)
        run_contrast(f"GSE220518_{key}_unadjusted", Wm, sm[col].values, [], foc); run_contrast(f"GSE220518_{key}_lineage_adjusted", Wm, sm[col].values, covariates(sm, lin), foc)
    m = mds & st.progression.isin(["MDS", "MDSp"]); Wm, sm = sub(m); y = (sm.progression == "MDSp").astype(float).values
    run_contrast("GSE220518_B1_progression_mutation_adjusted", Wm, y, covariates(sm, ["sf3b1", "srsf2", "u2af1"])); run_contrast("GSE220518_B1_progression_mutation_lineage_adjusted", Wm, y, covariates(sm, ["sf3b1", "srsf2", "u2af1"] + lin)); run_contrast("GSE220518_B1_progression_unadjusted", Wm, y, [])
    m2 = m & (st.nmut == 0); Wm, sm = sub(m2); run_contrast("GSE220518_B2_progression_within_SFwt_lineage_adjusted", Wm, (sm.progression == "MDSp").astype(float).values, covariates(sm, lin))
    run_contrast("GSE220518_D_MDS_vs_healthy_descriptive", W, mds.astype(float).values, [], ("SF3B1",))
elif gse == "GSE193946":
    S = marker_scores(list(st.index), EPI); st = st.join(S); cc = ["epithelial", "stroma", "immune"]; lev = st.histopathology.map({"BE.LGD": 0, "BE.HGD": 1, "EAC": 2}); print(dict(st.histopathology.value_counts())); print(st.groupby("histopathology")[list(S.columns)].mean().round(2).to_string())
    for key, m, y in (("HGD_vs_LGD", lev.isin([0, 1]), None), ("EAC_vs_LGD", lev.isin([0, 2]), None), ("stage_trend", lev.notna(), "trend")):
        Wm, sm = sub(m); grp = lev[m.values].values.astype(float) if y == "trend" else (lev[m.values].values > 0).astype(float)
        run_contrast(f"GSE193946_{key}_unadjusted", Wm, grp, []); run_contrast(f"GSE193946_{key}_composition_adjusted", Wm, grp, covariates(sm, cc))
elif gse == "GSE79209":
    S = marker_scores(list(st.index), EPI); st = st.join(S); cc = ["epithelial", "squamous", "ciliated", "immune"]; base = ["smoking_status", "sex", "copd_status"]; print(dict(st.dysplasia_status.value_counts(dropna=False)), dict(st.max_histology.value_counts())); print(st.groupby("dysplasia_status")[list(S.columns)].mean().round(2).to_string())
    m = st.dysplasia_status.isin(["Dysplasia", "Normal"]); Wm, sm = sub(m); y = (sm.dysplasia_status == "Dysplasia").astype(float).values
    run_contrast("GSE79209_dysplasia_vs_normal_adjusted", Wm, y, covariates(sm, base)); run_contrast("GSE79209_dysplasia_vs_normal_unadjusted", Wm, y, []); run_contrast("GSE79209_dysplasia_vs_normal_adjusted_plus_composition", Wm, y, covariates(sm, base + cc))
    lev = st.max_histology.map({"Normal": 0, "Hyperplasia": 1, "Metaplasia": 2, "Mild dysplasia": 3, "Moderate dysplasia": 4, "Severe dysplasia": 5}); m = lev.notna(); Wm, sm = sub(m)
    run_contrast("GSE79209_histology_trend_adjusted", Wm, lev[m.values].values.astype(float), covariates(sm, base)); run_contrast("GSE79209_histology_trend_adjusted_plus_composition", Wm, lev[m.values].values.astype(float), covariates(sm, base + cc))
elif gse == "GSE227919":
    S = marker_scores(list(st.index), EPI); st = st.join(S); cc = ["squamous", "stroma", "immune"]; lev = st.group.map({"Control": 0, "premalignant lesions": 1, "OSCC": 2}); print(dict(st.group.value_counts())); print(st.groupby("group")[list(S.columns)].mean().round(2).to_string())
    for key, m, y in (("premalignant_vs_control", lev.isin([0, 1]), None), ("OSCC_vs_control", lev.isin([0, 2]), None), ("stage_trend", lev.notna(), "trend")):
        Wm, sm = sub(m); grp = lev[m.values].values.astype(float) if y == "trend" else (lev[m.values].values > 0).astype(float)
        run_contrast(f"GSE227919_{key}_unadjusted", Wm, grp, []); run_contrast(f"GSE227919_{key}_composition_adjusted", Wm, grp, covariates(sm, cc))
elif gse == "GSE202048":
    S = marker_scores(list(st.index), EPI); st = st.join(S); cc = ["squamous", "stroma", "immune"]; tt = st.title.str.extract(r"^(\w+)\s+(\d+)\s*([MT])$"); st["kind"] = tt[0]; st["pid"] = tt[0] + "_" + tt[1]; st["side"] = tt[2]; print(st[["title", "pathology", "pid", "side"]].head(6).to_string()); print(dict(st.pathology.value_counts()))
    pairs = st.reset_index().pivot_table(index="pid", columns="side", values="gsm", aggfunc="first").dropna(); pairs["kind"] = pairs.index.str.split("_").str[0]; print("pairs", len(pairs), dict(pairs.kind.value_counts()))
    for key, pm in (("OPMD_lesion_vs_margin", pairs.kind == "OPMD"), ("SCC_lesion_vs_margin", pairs.kind == "SCC"), ("all_lesions_vs_margin", pairs.kind.notna())):
        pp = pairs[pm.values]; dc = st.loc[pp["T"], cc].values - st.loc[pp["M"], cc].values
        run_paired(f"GSE202048_{key}_paired_unadjusted", W, list(pp["T"]), list(pp["M"])); run_paired(f"GSE202048_{key}_paired_composition_adjusted", W, list(pp["T"]), list(pp["M"]), dcov=dc)
pd.DataFrame(summ).to_csv(f"{O}/h4b2_summary_{gse}.tsv", sep="\t", index=False, float_format="%.4f"); print("\n== summary =="); print(pd.DataFrame(summ).to_string(index=False)); print("DONE")
