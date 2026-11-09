"""Cell-composition marker scores per sample from STAR ReadsPerGene (reverse-stranded column), and a confounding check
of the significant rMATS events: group difference in scores, PSI ~ epithelial score correlation, and ΔPSI re-estimated
with epithelial score (+ batch) as covariates. Usage: python cell_composition_scores.py <cohort_dir> <rmats_dir>  (b1 of the rMATS run = case group)
Writes <cohort_dir>/cell_composition_scores.tsv and <rmats_dir>/summary_sig_events_composition_adjusted.tsv."""
import sys, re, numpy as np, pandas as pd, statsmodels.api as sm
cohort, rd = sys.argv[1:3]
GTF = "/well/blagden/users/eic789/LncRNASPLICE/gencode.v47.primary_assembly.annotation.gtf"
MARK = {"epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"], "ft_secretory": ["PAX8", "OVGP1", "KRT7"], "ft_ciliated": ["FOXJ1", "CAPS", "TPPP3", "PIFO"],
        "breast_luminal": ["KRT8", "KRT18", "ESR1", "GATA3"], "breast_basal": ["KRT5", "KRT14", "TP63"], "stroma": ["COL1A1", "COL1A2", "DCN", "LUM"],
        "smooth_muscle": ["ACTA2", "MYH11", "DES", "CNN1"], "immune": ["PTPRC"], "endothelial": ["PECAM1", "VWF", "CDH5"], "adipose": ["ADIPOQ", "LEP", "PLIN1"]}
st = pd.read_csv(f"{cohort}/sample_table.tsv", sep="\t").set_index("gsm")
sym = {}
for l in open(GTF):
    if l.startswith("#"): continue
    f = l.split("\t")
    if f[2] != "gene": continue
    sym[re.search('gene_id "([^"]+)"', f[8]).group(1)] = re.search('gene_name "([^"]+)"', f[8]).group(1)
def strand_col(path):
    # STAR ReadsPerGene columns: 1 unstranded, 2 forward, 3 reverse. A stranded library has far fewer N_noFeature reads in its own column.
    nf = pd.read_csv(path, sep="\t", header=None, nrows=4, index_col=0).loc["N_noFeature"]
    return 2 if nf[2] < nf[3] / 3 else 3 if nf[3] < nf[2] / 3 else 1
col = strand_col(f"{cohort}/star/{st.index[0]}.ReadsPerGene.out.tab"); print("count column", col, {1: "unstranded", 2: "forward", 3: "reverse"}[col])
cnt = {}
for g in st.index:
    t = pd.read_csv(f"{cohort}/star/{g}.ReadsPerGene.out.tab", sep="\t", header=None, index_col=0, skiprows=4)
    cnt[g] = t[col]
C = pd.DataFrame(cnt); C.index = C.index.map(sym); C = C.groupby(level=0).sum()
L = np.log2(C / C.sum() * 1e6 + 1)
S = pd.DataFrame({k: L.reindex(v).dropna().mean() for k, v in MARK.items()}).T
S = S.loc[[k for k in S.index if S.loc[k].notna().all()]]
scores = S.T; scores.index.name = "gsm"
scores = scores.join(st[["batch"]] if "batch" in st.columns else st[[]])
scores.to_csv(f"{cohort}/cell_composition_scores.tsv", sep="\t")
from scipy.stats import mannwhitneyu
gs = lambda p: re.search(r"(GSM\d+)", p).group(1)
b1 = [gs(p) for p in open(f"{rd}/b1.txt").read().strip().split(",")]; b2 = [gs(p) for p in open(f"{rd}/b2.txt").read().strip().split(",")]; samples = b1 + b2
scores["group"] = ["b1_case" if g in b1 else "b2_control" if g in b2 else "unused" for g in scores.index]
scores.to_csv(f"{cohort}/cell_composition_scores.tsv", sep="\t")
print(f"== {cohort}: marker scores (mean log2 CPM), b1 (n={len(b1)}) vs b2 (n={len(b2)})")
for k in S.index:
    a = scores.loc[b1, k]; b = scores.loc[b2, k]
    print(f"{k:15s} case {a.mean():5.2f}  rest {b.mean():5.2f}  diff {a.mean()-b.mean():+5.2f}  MWU p={mannwhitneyu(a, b).pvalue:.3g}")
# significant events: PSI vs epithelial score, covariate-adjusted dPSI
grp = np.array([1] * len(b1) + [0] * len(b2))
epi = scores.loc[samples, "epithelial"].values; batch = pd.get_dummies(scores.loc[samples, "batch"], drop_first=True).astype(float).values if "batch" in scores else np.zeros((len(samples), 0))
sig = pd.read_csv(f"{rd}/summary_sig_events.tsv", sep="\t"); out = []
for t in sig.type.unique():
    x = pd.read_csv(f"{rd}/{t}.MATS.JC.txt", sep="\t").set_index("ID")
    for i in sig[sig.type == t].ID:
        r = x.loc[i]; psi = np.array([float(v) if v != "NA" else np.nan for v in (r.IncLevel1 + "," + r.IncLevel2).split(",")])
        ok = ~np.isnan(psi)
        X = np.column_stack([np.ones(ok.sum()), grp[ok], epi[ok], batch[ok]])
        fit = sm.OLS(psi[ok], X).fit(); r_epi = np.corrcoef(psi[ok], epi[ok])[0, 1]
        out.append((t, i, r.geneSymbol, round(float(r.IncLevelDifference), 3), round(fit.params[1], 3), round(fit.pvalues[1], 4), round(r_epi, 3), round(fit.params[2], 3)))
o = pd.DataFrame(out, columns=["type", "ID", "gene", "dPSI_rmats", "dPSI_adj_epi_batch", "p_group_adj", "r_psi_epithelial", "beta_epithelial"])
o.to_csv(f"{rd}/summary_sig_events_composition_adjusted.tsv", sep="\t", index=False)
print(f"\n== {len(o)} significant events: group effect after adjusting for epithelial score{' + batch' if batch.shape[1] else ''}")
print("sign kept:", (np.sign(o.dPSI_rmats) == np.sign(o.dPSI_adj_epi_batch)).mean().round(2), "| |adj dPSI|>0.1:", (o.dPSI_adj_epi_batch.abs() > 0.1).mean().round(2), "| p_group_adj<0.05:", (o.p_group_adj < 0.05).mean().round(2))
print("median |r(PSI, epithelial)|:", o.r_psi_epithelial.abs().median().round(2), "| events with |r|>0.5:", (o.r_psi_epithelial.abs() > 0.5).sum())
print(o.sort_values("p_group_adj").head(12).to_string(index=False))
