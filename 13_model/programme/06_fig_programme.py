"""Figure for the programme-level model: (a) subspace R2 by descriptor set, (b, c) programme-1 agreement by three metrics, (d) programmes passing.
Reads results_summary_final.tsv and results_programme_bootstrap_final.tsv. Output: 13_model/programme/fig_programme.png/.pdf"""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib import font_manager
for f in ("Arial.ttf", "Arial-Bold.ttf"): font_manager.fontManager.addfont(f"/well/blagden/users/eic789/.fonts/{f}")
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "savefig.dpi": 300, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
O = "/well/blagden/users/eic789/PerturbBridge/13_model/programme"; S = pd.read_csv(f"{O}/results_summary_final.tsv", sep="\t"); B = pd.read_csv(f"{O}/results_programme_bootstrap_final.tsv", sep="\t")
col = {"HepG2": "#2f5d8a", "K562": "#8fb4d9"}; cells = ("HepG2", "K562"); fig, ax = plt.subplots(2, 2, figsize=(7.09, 5.0)); ax = ax.ravel()
def get(c, s, m="ridge"): r = S[(S.cell == c) & (S.set == s) & (S.method == m)]; return r.iloc[0]
sets = [("S0", "sequence\n+ text"), ("ps_only", "Perturb-seq\nalone"), ("dep_only", "DepMap\nalone"), ("net_only", "network\nalone"), ("seq_net", "seq + text\n+ network\n(selected)"), ("seq_net_dep", "seq + text\n+ network\n+ DepMap")]
for i, c in enumerate(cells):
    x = np.arange(len(sets)) + (i - 0.5) * 0.36; v = [get(c, s).r2_sub for s, _ in sets]; ax[0].bar(x, v, 0.34, color=col[c], label=c)
    ax[0].scatter(x, [get(c, s).r2_sub_perm_p95 for s, _ in sets], marker="_", s=40, color="k", linewidths=0.8, zorder=3, label="permuted descriptors (95th pct)" if i == 0 else None)
    ax[0].axhline(get(c, "oracle_other_cell", "-").r2_sub, color=col[c], ls="--", lw=0.8)
ax[0].text(len(sets) - 0.45, get("K562", "oracle_other_cell", "-").r2_sub + 0.004, "same RBP measured in the\nother cell line (dashed)", ha="right", va="bottom", fontsize=5.5)
ax[0].set_xticks(range(len(sets))); ax[0].set_xticklabels([l for _, l in sets], fontsize=5.8); ax[0].set_ylabel("Variance explained in the 10-programme\nsubspace, held-out RBPs (R²)"); ax[0].legend(frameon=False, fontsize=5.8, loc="upper left"); ax[0].set_title("a  Descriptor sets (ridge)", loc="left", fontweight="bold", fontsize=7.5)
mets = [("r1", "Pearson"), ("rho1", "Spearman"), ("r1_trim", "Pearson, top 5 %\nscores removed")]; preds = [("S0", "ridge", "sequence + text", "#b8b8b8"), ("seq_net", "ridge", "+ network (selected)", "#3a8f5c"), ("seq_net_dep", "ridge", "+ network + DepMap", "#9fd0b1"), ("oracle_other_cell", "-", "measured in other cell line", "#6b4c9a")]
for i, c in enumerate(cells):
    a = ax[1 + i]
    for j, (s, m, lab, cc) in enumerate(preds): a.bar(np.arange(3) + (j - 1.5) * 0.2, [get(c, s, m)[k] for k, _ in mets], 0.19, color=cc, label=lab)
    a.set_xticks(range(3)); a.set_xticklabels([l for _, l in mets], fontsize=6); a.set_ylim(0, 0.72); a.set_ylabel("Programme-1 score, predicted vs measured"); a.set_title(f"{'bc'[i]}  Programme 1 (core spliceosome axis), {c}", loc="left", fontweight="bold", fontsize=7.5)
    if i == 0: a.legend(frameon=False, fontsize=5.8, loc="upper right")
ps = [("S0", "sequence\n+ text"), ("net_only", "network\nalone"), ("seq_net", "sequence + text\n+ network"), ("seq_net_dep", "+ DepMap")]
for i, c in enumerate(cells):
    n = [int(B[(B.cell == c) & (B.set == s) & (B.method == "ridge") & (B.programme > 1)].passed.sum()) for s, _ in ps]; x = np.arange(len(ps)) + (i - 0.5) * 0.36; ax[3].bar(x, n, 0.34, color=col[c], label=c)
    for xx, nn in zip(x, n): ax[3].text(xx, nn + 0.15, str(nn), ha="center", fontsize=6)
ax[3].set_xticks(range(len(ps))); ax[3].set_xticklabels([l for _, l in ps], fontsize=5.8); ax[3].set_ylim(0, 10); ax[3].set_ylabel("Programmes 2–10 predicted above\npermuted controls (of 9)"); ax[3].legend(frameon=False, fontsize=5.8); ax[3].set_title("d  Programmes beyond the first", loc="left", fontweight="bold", fontsize=7.5)
fig.tight_layout(); fig.savefig(f"{O}/fig_programme.png"); fig.savefig(f"{O}/fig_programme.pdf"); print("saved")
