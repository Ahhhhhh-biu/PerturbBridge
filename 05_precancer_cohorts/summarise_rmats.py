"""Summarise a group-level rMATS output directory: coverage-filtered event counts, significant events, top SE genes.
Usage: python summarise_rmats.py <rmats_dir> <label>. Writes <rmats_dir>/summary_sig_events.tsv."""
import sys, numpy as np, pandas as pd
d, label = sys.argv[1], sys.argv[2]
def cov(col): return col.str.split(",").apply(lambda v: np.mean([float(x) for x in v if x != "NA"]) if any(x != "NA" for x in v) else 0)
rows, sig_all = [], []
for t in ["SE", "A3SS", "A5SS", "MXE", "RI"]:
    x = pd.read_csv(f"{d}/{t}.MATS.JC.txt", sep="\t")
    x["cov1"] = cov(x.IJC_SAMPLE_1) + cov(x.SJC_SAMPLE_1); x["cov2"] = cov(x.IJC_SAMPLE_2) + cov(x.SJC_SAMPLE_2)
    x = x[(x.cov1 >= 10) & (x.cov2 >= 10)]
    s = x[(x.FDR < 0.05) & (x.IncLevelDifference.abs() > 0.1)].copy(); s["type"] = t
    rows.append((t, len(x), len(s), s.geneSymbol.nunique()))
    sig_all.append(s[["type", "ID", "GeneID", "geneSymbol", "chr", "strand", "PValue", "FDR", "IncLevelDifference", "cov1", "cov2"]])
sig = pd.concat(sig_all).sort_values("FDR"); sig.to_csv(f"{d}/summary_sig_events.tsv", sep="\t", index=False)
print(f"== {label}: events with mean(IJC+SJC)>=10 in both groups; significant = FDR<0.05 & |dPSI|>0.1")
print(pd.DataFrame(rows, columns=["type", "n_cov10", "n_sig", "n_sig_genes"]).to_string(index=False))
print("total sig:", len(sig), "| direction: dPSI>0 (b1 higher)", (sig.IncLevelDifference > 0).sum(), "dPSI<0", (sig.IncLevelDifference < 0).sum())
print(sig[sig.type == "SE"].head(15)[["geneSymbol", "FDR", "IncLevelDifference", "cov1", "cov2"]].round(4).to_string(index=False))
