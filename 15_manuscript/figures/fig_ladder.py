"""Ladder figure: for each task the chance level or control, the model, and the ceiling or reference, on the task's own scale.
Sources: event level 13_model/benchmark/master_table.tsv (v3 deep, v3 shuffled RBP) and 13_model/clean_ceiling_cross_run.tsv (same-RBP independent
experiments, mean sign_agree); programme level 13_model/programme/results_summary_final.tsv (set seq_net, ridge; permuted 95th percentile; oracle =
same RBP measured in the other cell line); ENCODE retrieval 14_inverse/h3_encode_summary.tsv (setting cross); TCGA 14_inverse/h3_tcga_tests.tsv
(tests with the mutated gene in the dictionary); premalignant cohorts 14_inverse/h4_summary.tsv and h4new_summary.tsv (adjusted analyses)."""
import numpy as np, pandas as pd
from _style import *
from matplotlib.lines import Line2D
M = pd.read_csv(f"{PB}/13_model/benchmark/master_table.tsv", sep="\t").set_index(["model", "cell_line"]).sign_acc; C = pd.read_csv(f"{PB}/13_model/clean_ceiling_cross_run.tsv", sep="\t").groupby("cell_line").sign_agree.mean()
S = pd.read_csv(f"{PB}/13_model/programme/results_summary_final.tsv", sep="\t").set_index(["cell", "set", "method"]); E = pd.read_csv(f"{PB}/14_inverse/h3_encode_summary.tsv", sep="\t"); E = E[E.setting == "cross"].set_index(["query_cell", "repro_class"])
T = pd.read_csv(f"{PB}/14_inverse/h3_tcga_tests.tsv", sep="\t"); T = T[T.gene_rank.notna()].copy(); T["frac"] = T.gene_rank / T.n_rbps
H = pd.read_csv(f"{PB}/14_inverse/h4_summary.tsv", sep="\t"); H = H[H.adjustment == "adjusted"]; N = pd.read_csv(f"{PB}/14_inverse/h4new_summary.tsv", sep="\t"); N = N[N.contrast.str.startswith("GSE182338") & N.contrast.str.endswith("_adjusted")]
rows = [("Single events, unseen RBP", "direction accuracy on clean events", (0.47, 0.60), {c: (M[("v3 shuffled RBP", c)], M[("v3 deep", c)], C[c]) for c in CELLS}, 0.5),
        ("Programme 1, unseen RBP", "Pearson r of predicted vs measured score", (-0.02, 0.70), {c: (S.loc[(c, "seq_net", "ridge"), "r1_perm_p95"], S.loc[(c, "seq_net", "ridge"), "r1"], S.loc[(c, "oracle_other_cell", "-"), "r1"]) for c in CELLS}, None),
        ("Ten programmes, unseen RBP", "variance explained (R²)", (-0.01, 0.20), {c: (S.loc[(c, "seq_net", "ridge"), "r2_sub_perm_p95"], S.loc[(c, "seq_net", "ridge"), "r2_sub"], S.loc[(c, "oracle_other_cell", "-"), "r2_sub"]) for c in CELLS}, None),
        ("Inverse mapping, ENCODE", "top-5 hit rate across cell lines", (-0.03, 0.90), {c: (E.loc[(c, "all"), "expected_top5"], E.loc[(c, "all"), "top5"], E.loc[(c, "high"), "top5"]) for c in CELLS}, None)]
fig = plt.figure(figsize=(7.09, 4.5)); n = len(rows) + 2; L, R, top, bot = 0.36, 0.975, 0.92, 0.03; hgt = (top - bot) / n
for i, (name, sub, xl, vals, chance) in enumerate(rows):
    ax = fig.add_axes([L, top - (i + 1) * hgt + 0.045, R - L, hgt - 0.065]); ax.set_ylim(-0.6, 0.6); ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_xlim(*xl)
    for j, c in enumerate(CELLS):
        y = 0.25 - 0.5 * j; ctl, mod, ceil = vals[c]; ax.plot([min(ctl, mod, ceil), max(ctl, mod, ceil)], [y, y], color="0.8", lw=0.8, zorder=1)
        ax.scatter([ctl], [y], s=24, facecolor="white", edgecolor="0.4", lw=0.8, zorder=2); ax.scatter([ceil], [y], s=26, marker="D", color="k", zorder=2, lw=0); ax.scatter([mod], [y], s=30, color=COL[c], zorder=3, edgecolor="white", lw=0.4)
    if chance is not None: ax.axvline(chance, color="k", lw=0.5, ls="--")
    fig.text(0.015, top - (i + 0.5) * hgt + 0.012, name, fontsize=7.5, fontweight="bold", va="center"); fig.text(0.015, top - (i + 0.5) * hgt - 0.016, sub, fontsize=6, color="0.3", va="center")
fig.text(0.015 + 0.0, top - 3.5 * hgt - 0.042, "reference: high-reproducibility RBPs", fontsize=5.5, color="0.4", va="center")
fig.text(0.015, top - 1.5 * hgt - 0.042, "reference: same RBP measured in the other cell line", fontsize=5.5, color="0.4", va="center"); fig.text(0.015, top - 2.5 * hgt - 0.042, "reference: same RBP measured in the other cell line", fontsize=5.5, color="0.4", va="center"); fig.text(0.015, top - 0.5 * hgt - 0.042, "ceiling: same RBP, independent experiments", fontsize=5.5, color="0.4", va="center")
# TCGA
i = len(rows); ax = fig.add_axes([L, top - (i + 1) * hgt + 0.045, R - L, hgt - 0.065]); rng = np.random.default_rng(0); ax.scatter(T.frac, rng.uniform(-0.3, 0.3, len(T)), s=8, color="0.6", lw=0, zorder=2); ax.axvline(0.5, color="k", lw=0.5, ls="--")
ax.scatter([T.frac.median()], [0], s=34, color="#6b4c9a", zorder=4, edgecolor="white", lw=0.4); u = T[T.test == "TCGA-UVM_SF3B1_hotspot"].iloc[0]; ax.scatter([u.frac], [0.0], s=20, color="#c0392b", zorder=5, lw=0); ax.annotate("SF3B1 hotspot, uveal melanoma", (u.frac, 0.0), (0.06, 0.42), fontsize=5.6, color="#c0392b", va="center", arrowprops=dict(arrowstyle="-", lw=0.4, color="#c0392b"))
ax.text(T.frac.median(), -0.52, f"median {T.frac.median():.2f}", fontsize=5.6, color="#6b4c9a", ha="right", va="center"); ax.set_xlim(-0.03, 1.0); ax.set_ylim(-0.7, 0.6); ax.set_yticks([]); ax.spines["left"].set_visible(False)
fig.text(0.015, top - (i + 0.5) * hgt + 0.012, "Inverse mapping, TCGA tumours", fontsize=7.5, fontweight="bold", va="center"); fig.text(0.015, top - (i + 0.5) * hgt - 0.016, f"rank fraction of the mutated RBP, {len(T)} tests", fontsize=6, color="0.3", va="center"); fig.text(0.015, top - (i + 0.5) * hgt - 0.042, "lower is better; 0.5 expected by chance", fontsize=5.5, color="0.4", va="center")
# premalignant
i += 1; ax = fig.add_axes([L, top - (i + 1) * hgt + 0.045, R - L, hgt - 0.065]); pv = [(r.p_top_vs_maxnull, "#b5483a") for r in H.itertuples()] + [(r.p_top_vs_maxnull, "#d9a49c") for r in N.itertuples()]
ax.scatter([p for p, _ in pv], rng.uniform(-0.25, 0.25, len(pv)), s=22, color=[c for _, c in pv], zorder=3, edgecolor="white", lw=0.4); ax.axvline(0.05, color="k", lw=0.5, ls="--"); ax.text(0.058, 0.47, "0.05", fontsize=5.6, va="center"); ax.set_xlim(-0.03, 1.0); ax.set_ylim(-0.6, 0.6); ax.set_yticks([]); ax.spines["left"].set_visible(False)
fig.text(0.015, top - (i + 0.5) * hgt + 0.012, "Inverse mapping, premalignant tissue", fontsize=7.5, fontweight="bold", va="center"); fig.text(0.015, top - (i + 0.5) * hgt - 0.016, "p of the top RBP against the permutation maximum", fontsize=6, color="0.3", va="center")
fig.text(0.015, top - (i + 0.5) * hgt - 0.052, f"adjusted analyses: fallopian tube and breast (dark);\nsorted mammary cells, {len(N)} contrasts (light)", fontsize=5.5, color="0.4", va="center", linespacing=1.2)
hd = [Line2D([], [], marker="o", ls="", mfc="white", mec="0.4", ms=5, label="chance or control"), Line2D([], [], marker="o", ls="", color=COL["HepG2"], ms=5.5, label="model, HepG2"), Line2D([], [], marker="o", ls="", color=COL["K562"], ms=5.5, label="model, K562"), Line2D([], [], marker="D", ls="", color="k", ms=4.5, label="ceiling or reference")]
fig.legend(handles=hd, loc="upper right", bbox_to_anchor=(R, 0.995), ncol=4, fontsize=6.3, handletextpad=0.2, columnspacing=1.2)
save(fig, "fig_ladder")
