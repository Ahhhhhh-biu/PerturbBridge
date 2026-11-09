"""Block 2b: H3 on public tumours. For every (TCGA project, RBP) with >= MINMUT mutant samples that have Snaptron junction data, compare
per-sample PSI (14_inverse/tcga_psi/part_*.parquet, coverage >= 10, events covered in >= 70% of the project's samples) between mutant
and wild-type primary tumours (Welch t), take the query (p < 0.01 & |dPSI| > 0.05, capped 500) and rank the mutated RBP in the
dictionary. Mutation classes: truncating (nonsense / frameshift / splice site) for all RBPs; hotspot missense for SF3B1 (R625*, K666*,
K700E), U2AF1 (S34*, Q157*), SRSF2 (P95*). RBM10 (31 truncating LUAD) is not an ENCODE target: reported as 'what does it retrieve'.
Null: 200 permutations of the mutant label within the project. Also tumour vs matched normal (BRCA, LUAD) as a descriptive query.
Outputs: 14_inverse/h3_tcga_tests.tsv (one row per test), h3_tcga_rankings.tsv (top 15 per test); stdout -> logs/h3_tcga.txt"""
import sys, glob, numpy as np, pandas as pd; from scipy import stats; sys.path.insert(0, "/well/blagden/users/eic789/PerturbBridge/14_inverse"); import retrieval as R
PB = R.PB; O = R.O; rng = np.random.default_rng(0); pd.set_option("display.width", 250); MINMUT = 5; NPERM = 200
P = R.load_dictionary(); info = R.rbp_info()
S = pd.read_csv(f"{O}/tcga_samples.tsv", sep="\t", dtype=str, usecols=["rail_id", "tcga_barcode", "gdc_cases.project.project_id", "cgc_sample_sample_type"]).rename(columns={"gdc_cases.project.project_id": "project", "cgc_sample_sample_type": "stype"})
S["rail_id"] = S.rail_id.astype(int); S["sid"] = S.tcga_barcode.str[:15]; S = S.drop_duplicates("rail_id")
mut = pd.read_csv(f"{O}/tcga_mutations.tsv", sep="\t"); seq = pd.read_csv(f"{O}/tcga_sequenced_samples.tsv", sep="\t"); study2proj = {"uvm": "TCGA-UVM", "luad": "TCGA-LUAD", "brca": "TCGA-BRCA", "ov": "TCGA-OV", "laml": "TCGA-LAML", "skcm": "TCGA-SKCM"}
mut["project"] = mut.study.str.split("_").str[0].map(study2proj); seq["project"] = seq.study.str.split("_").str[0].map(study2proj)
TRUNC = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"}
HOT = {"SF3B1": r"^(R625|K666|K700|E622|H662|G742)", "U2AF1": r"^(S34|Q157)", "SRSF2": r"^P95"}
mut["cls"] = np.where(mut.mutationType.isin(TRUNC), "truncating", np.where([bool(__import__("re").match(HOT.get(g, r"$^"), str(pc))) for g, pc in zip(mut.gene, mut.proteinChange)], "hotspot", "other"))
parts = sorted(glob.glob(f"{O}/tcga_psi/part_*.parquet")); X = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True); print("psi rows", len(X), "events", X.event_id.nunique(), "parts", len(parts))
X = X[X.inc + X["skip"] >= 10].merge(S[["rail_id", "sid", "project", "stype"]], on="rail_id")
tests, ranks = [], []
def run_test(label, project, Wp, is_case, n_case_label):
    Y = Wp.values; a = Y[is_case]; b = Y[~is_case]
    def query_from(mask):
        A = Y[mask]; B = Y[~mask]; t, p = stats.ttest_ind(A, B, axis=0, equal_var=False, nan_policy="omit"); d = np.nanmean(A, 0) - np.nanmean(B, 0)
        sel = np.where((p < 0.01) & (np.abs(d) > 0.05))[0]
        if len(sel) > 500: sel = sel[np.argsort(p[sel])[:500]]
        return pd.Series(d[sel], index=Wp.columns[sel]), d, p
    q, d_all, p_all = query_from(is_case)
    if len(q) < 50: tests.append(dict(test=label, project=project, n_case=int(is_case.sum()), n_control=int((~is_case).sum()), n_query=len(q), note="too few query events")); print(label, "too few query events", len(q)); return
    s, n = R.score_against(q, P); rk = R.rank_report(s, info, label)
    null = np.full((NPERM, len(P.columns)), np.nan)
    for k in range(NPERM):
        qk, _, _ = query_from(rng.permutation(is_case))
        if len(qk) >= 50: sk, _ = R.score_against(qk, P); null[k] = sk.reindex(P.columns).values
    obs = s.reindex(P.columns).values; emp = np.array([(np.nansum(null[:, j] >= obs[j]) + 1) / (np.isfinite(null[:, j]).sum() + 1) if np.isfinite(obs[j]) else np.nan for j in range(len(obs))])
    rk = rk.merge(pd.DataFrame({"rbp": P.columns, "perm_p": emp, "null_mean": np.nanmean(null, 0), "null_sd": np.nanstd(null, 0)}), on="rbp"); rk["z_vs_null"] = (rk.score - rk.null_mean) / rk.null_sd
    gene = n_case_label; row = dict(test=label, project=project, n_case=int(is_case.sum()), n_control=int((~is_case).sum()), n_events_tested=int(np.isfinite(p_all).sum()), n_query=len(q), n_rbps=int(s.notna().sum()), top1=rk.rbp.iloc[0], top1_score=rk.score.iloc[0], top1_perm_p=rk.perm_p.iloc[0],
               top5=";".join(rk.rbp.iloc[:5]), n_rbp_perm_p05=int((rk.perm_p < 0.05).sum()))
    if gene in rk.rbp.values: r0 = rk[rk.rbp == gene].iloc[0]; row.update(gene=gene, gene_rank=int(r0["rank"]), gene_score=r0.score, gene_perm_p=r0.perm_p, gene_repro=r0.repro_class)
    else: row.update(gene=gene, gene_rank=np.nan, note="gene not in dictionary")
    tests.append(row); ranks.append(rk.head(15)); print(f"{label}: cases {is_case.sum()} query {len(q)} | gene rank {row.get('gene_rank')} score {row.get('gene_score', np.nan):.3f} perm_p {row.get('gene_perm_p', np.nan):.3f} | top5 {row['top5']}")
for project, Xp in X.groupby("project"):
    prim = Xp[Xp.stype.str.contains("Primary", na=False) | (project == "TCGA-SKCM")]   # SKCM is mostly metastatic
    W = prim.pivot_table(index="sid", columns="event_id", values="psi", aggfunc="mean"); W = W.loc[:, W.notna().mean(0) >= 0.7]
    sequenced = set(seq.loc[seq.project == project, "sampleId"]); W = W.loc[W.index.isin(sequenced)]; print(f"\n== {project}: primary samples with junctions and mutation calls {W.shape[0]}, events {W.shape[1]}")
    mp = mut[(mut.project == project) & mut.cls.isin(["truncating", "hotspot"])]
    for (gene, cls), g in mp.groupby(["gene", "cls"]):
        cases = set(g.sampleId) & set(W.index)
        if len(cases) < MINMUT: continue
        others = set(mut.loc[(mut.project == project) & (mut.gene == gene), "sampleId"])   # any mutation in the gene -> exclude from wild type
        keep = W.index.isin(cases) | ~W.index.isin(others); Wp = W.loc[keep]; is_case = Wp.index.isin(cases)
        run_test(f"{project}_{gene}_{cls}", project, Wp, is_case, gene)
    # tumour vs normal
    norm = Xp[Xp.stype.str.contains("Normal", na=False)]
    if norm.sid.nunique() >= 10:
        Wn = Xp[Xp.stype.str.contains("Primary|Normal", na=False)].pivot_table(index=["sid", "stype"], columns="event_id", values="psi", aggfunc="mean"); Wn = Wn.loc[:, Wn.notna().mean(0) >= 0.7]
        is_t = np.array([("Primary" in st) for _, st in Wn.index]); Wn.index = [f"{a}|{b}" for a, b in Wn.index]; run_test(f"{project}_tumour_vs_normal", project, Wn, is_t, "none")
T = pd.DataFrame(tests); T.to_csv(f"{O}/h3_tcga_tests.tsv", sep="\t", index=False, float_format="%.4f"); pd.concat(ranks).to_csv(f"{O}/h3_tcga_rankings.tsv", sep="\t", index=False, float_format="%.4f") if ranks else None
print("\n== summary ==\n", T.round(3).to_string(index=False)); print("DONE")
