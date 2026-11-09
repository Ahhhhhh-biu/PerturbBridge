"""Block 4b (post-freeze in-house validation): independent-layer consistency of the H4 fallopian-tube RBP set in the TMA221/222
spatial proteomics (epithelial ROIs; NM normal mucosa -> ST STIC precursor -> TM tumour). If knockdown of RBP X mimics the carrier
tissue splicing shift (positive dictionary score), X's protein is expected to DECREASE along the progression axis.
Set-level test only (single-protein claims were shown not to survive abundance matching in LARP1_Biomarker_paper): signed mean
log2FC (ST vs NM, TM vs NM) of the selected set vs 2,000 abundance- and detection-matched random protein sets of the same size.
Sets tested: H4 adjusted ranking perm_p < 0.05 (fallback: top 10, flagged); all high-reproducibility dictionary RBPs; ribosomal
proteins (positive control for sensitivity, direction-free). Outputs: 14_inverse/tma_set_consistency.tsv, tma_selected_proteins.tsv;
stdout -> logs/tma_consistency.txt"""
import numpy as np, pandas as pd
from scipy import stats
PB = "/well/blagden/users/eic789/PerturbBridge"; O = f"{PB}/14_inverse"; RAW = "/well/blagden/users/eic789/Preliminary spatial proteomics data of TMA221 and 222/raw_data"; rng = np.random.default_rng(7); pd.set_option("display.width", 250)
norm = pd.read_csv(f"{RAW}/TMA221_TMA222_SpatialProt_Epi_norm.csv", index_col=0); meta = pd.read_csv(f"{RAW}/TMA221_TMA222_SpatialProt_Epi_Metadata.csv").set_index("sample_id").reindex(norm.columns)
sym = pd.read_csv(f"{RAW}/TMA221_TMA222_SpatialProt_Epi_ProteinMapping.csv"); g2p = {}
for p, g in zip(sym["Protein.Group"], sym.Genes.astype(str)):
    if p not in norm.index: continue
    for x in str(g).split(";"): g2p.setdefault(x, p)
tt = meta.tissue_type; cols = {k: norm.columns[(tt == k).values] for k in ("NM", "ST", "TM")}; print("ROIs:", {k: len(v) for k, v in cols.items()}, "| proteins", norm.shape[0])
lfc = pd.DataFrame({"ST_vs_NM": norm[cols["ST"]].mean(1) - norm[cols["NM"]].mean(1), "TM_vs_NM": norm[cols["TM"]].mean(1) - norm[cols["NM"]].mean(1)}); ab = norm.mean(1); det = norm.notna().mean(1)
rk = pd.read_csv(f"{O}/h4_GSE274628_fallopian_adjusted_ranking.tsv", sep="\t"); info = pd.read_csv(f"{O}/dict_rbps.tsv", sep="\t")
sel = rk[rk.perm_p < 0.05]; note = "perm_p<0.05"
if len(sel) < 3: sel = rk.head(10); note = "fallback: top 10 by score (no RBP reached perm_p<0.05)"
sets = {"H4_fallopian_selected": (sel.rbp.tolist(), np.sign(sel.score.values)), "dictionary_high_repro": (sorted(set(info[info.repro_class == "high"].target)), None),
        "ribosomal_positive_control": ([str(g) for g in sym.Genes if str(g).startswith(("RPL", "RPS")) and str(g).isalnum()], None)}
def matched_null(prots, n=2000):
    pool = norm.index[~norm.index.isin(prots)]; A = ab[pool].values; D = det[pool].values; out = []
    for _ in range(n):
        pick = []
        for p in prots:
            cand = pool[(np.abs(A - ab[p]) < 0.5) & (np.abs(D - det[p]) < 0.1)]; pick.append(rng.choice(cand) if len(cand) else rng.choice(pool))
        out.append(pick)
    return out
rows, prot_rows = [], []
for name, (genes, signs) in sets.items():
    prots = [g2p[g] for g in genes if g in g2p]; ok = [i for i, g in enumerate(genes) if g in g2p]; sg = signs[ok] if signs is not None else np.ones(len(prots))
    if len(prots) < 3: print(name, "too few proteins measured", len(prots)); continue
    for g, p, s in zip([genes[i] for i in ok], prots, sg): prot_rows.append(dict(set=name, gene=g, protein=p, expected_sign=-s if signs is not None else np.nan, ST_vs_NM=lfc.ST_vs_NM[p], TM_vs_NM=lfc.TM_vs_NM[p], abundance=ab[p], detection=det[p]))
    null = matched_null(prots)
    for ax in ("ST_vs_NM", "TM_vs_NM"):
        v = lfc[ax][prots].values; obs = np.nanmean(-sg * v) if signs is not None else np.nanmean(np.abs(v))   # signed: expected decrease -> positive; unsigned sets: magnitude
        nv = np.array([np.nanmean(-sg * lfc[ax][pk].values) if signs is not None else np.nanmean(np.abs(lfc[ax][pk].values)) for pk in null])
        rows.append(dict(set=name, n_genes=len(genes), n_measured=len(prots), axis=ax, statistic="signed mean log2FC (expected direction positive)" if signs is not None else "mean |log2FC|", observed=obs, null_mean=nv.mean(), null_p95=np.percentile(nv, 95), empirical_p=(np.sum(nv >= obs) + 1) / (len(nv) + 1),
                         frac_expected_direction=float(np.mean(np.sign(v) == -sg)) if signs is not None else np.nan, note=note if name.startswith("H4") else ""))
T = pd.DataFrame(rows); T.to_csv(f"{O}/tma_set_consistency.tsv", sep="\t", index=False, float_format="%.4f"); pd.DataFrame(prot_rows).to_csv(f"{O}/tma_selected_proteins.tsv", sep="\t", index=False, float_format="%.4f")
print(T.round(3).to_string(index=False)); print("\nselected set:", note); print(pd.DataFrame(prot_rows).query("set == 'H4_fallopian_selected'").round(3).to_string(index=False)); print("DONE")
