"""Figure: junction-level (rMATS) shRNA-vs-CRISPR ΔPSI agreement per RBP vs TPM/SUPPA-level and gene-level Δlog2TPM."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
D = "/well/blagden/users/eic789/PerturbBridge"
j = pd.read_csv(f"{D}/08_junction_psi/junction_cross_technology.tsv", sep="\t")
t = pd.read_csv(f"{D}/04_h2_check/h2b_cross_technology.tsv", sep="\t")
tp = t.set_index("target")[["dPSI_cross_tech_corr", "dlog2TPM_cross_tech_corr"]].rename(columns={"dPSI_cross_tech_corr": "dPSI", "dlog2TPM_cross_tech_corr": "dlog2TPM"})
j = j.sort_values("r_cov50", ascending=False).reset_index(drop=True)
x = np.arange(len(j)); w = 0.22
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(7.2, 3.2))
ax.bar(x - 1.5*w, j.r_all, w, color="#9ecae1", label="junction ΔPSI, all events")
ax.bar(x - 0.5*w, j.r_cov50, w, color="#3182bd", label="junction ΔPSI, coverage ≥ 50")
ax.bar(x + 0.5*w, j.r_sig_either, w, color="#08519c", label="junction ΔPSI, significant in either")
if True:
    ax.bar(x + 1.5*w, [tp.loc[r, "dPSI"] if r in tp.index else np.nan for r in j.target], w, color="#fdae6b", label="TPM-derived ΔPSI (SUPPA2)")
    ax.plot(x, [tp.loc[r, "dlog2TPM"] if r in tp.index else np.nan for r in j.target], "k_", ms=10, mew=1.5, label="gene-level Δlog2TPM")
ax.plot(x, j.null_median_r, "r_", ms=10, mew=1.2, label="random-pair null (median)")
ax.set_xticks(x); ax.set_xticklabels([f"{r.target}\n{r.cell_line}" for r in j.itertuples()])
ax.set_ylabel("shRNA vs CRISPR Pearson r"); ax.axhline(0, color="grey", lw=0.5)
ax.set_title("Same-RBP cross-technology agreement, ENCODE knockdowns", fontsize=9)
ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper right")
fig.tight_layout(); fig.savefig(f"{D}/08_junction_psi/junction_cross_technology.png", dpi=300); fig.savefig(f"{D}/08_junction_psi/junction_cross_technology.pdf")
print("saved")
