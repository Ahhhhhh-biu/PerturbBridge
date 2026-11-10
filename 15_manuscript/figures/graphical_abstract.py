"""Graphical abstract (180 x 90 mm): compendium with two artefacts removed -> forward prediction for unseen RBPs -> inverse mapping.
Numbers are read from repository tables: runs 13_model/ceiling_within_run.tsv, RBPs 11_splits/rbp_folds.tsv, control-deviation rule
13_model/clean_target_summary.tsv, event level 13_model/benchmark/master_table.tsv and 13_model/clean_ceiling_cross_run.tsv, programme level
13_model/programme/results_summary_final.tsv, retrieval 14_inverse/h3_encode_summary.tsv, tumours 14_inverse/h3_tcga_tests.tsv, premalignant
cohorts 14_inverse/h4_summary.tsv and h4new_summary.tsv (adjusted analyses)."""
import numpy as np, pandas as pd
from _style import *
from matplotlib.patches import FancyBboxPatch, Ellipse, Rectangle, Circle, FancyArrowPatch
n_runs = len(pd.read_csv(f"{PB}/13_model/ceiling_within_run.tsv", sep="\t")); n_rbp = len(pd.read_csv(f"{PB}/11_splits/rbp_folds.tsv", sep="\t"))
CT = pd.read_csv(f"{PB}/13_model/clean_target_summary.tsv", sep="\t").set_index(["cell_line", "set"]).rule; rule_sig = CT.xs("sig_FDR", level=1).mean(); rule_clean = CT.xs("clean_0.05", level=1).mean()
M = pd.read_csv(f"{PB}/13_model/benchmark/master_table.tsv", sep="\t").set_index(["model", "cell_line"]).sign_acc; ev = [M[("v3 deep", c)] for c in CELLS]; ceil = pd.read_csv(f"{PB}/13_model/clean_ceiling_cross_run.tsv", sep="\t").groupby("cell_line").sign_agree.mean()
S = pd.read_csv(f"{PB}/13_model/programme/results_summary_final.tsv", sep="\t").set_index(["cell", "set", "method"]); r1 = [S.loc[(c, "seq_net", "ridge"), "r1"] for c in CELLS]; ref = np.mean([S.loc[(c, "oracle_other_cell", "-"), "r1"] for c in CELLS]); perm = np.mean([S.loc[(c, "seq_net", "ridge"), "r1_perm_p95"] for c in CELLS])
E = pd.read_csv(f"{PB}/14_inverse/h3_encode_summary.tsv", sep="\t"); E = E[(E.setting == "cross") & (E.repro_class == "all")]; top5 = E.top5.mean(); exp5 = E.expected_top5.mean()
T = pd.read_csv(f"{PB}/14_inverse/h3_tcga_tests.tsv", sep="\t"); T = T[T.gene_rank.notna()]; H = pd.read_csv(f"{PB}/14_inverse/h4_summary.tsv", sep="\t"); N = pd.read_csv(f"{PB}/14_inverse/h4new_summary.tsv", sep="\t"); N = N[N.contrast.str.startswith("GSE182338") & N.contrast.str.endswith("_adjusted")]
padj = list(H[H.adjustment == "adjusted"].p_top_vs_maxnull) + list(N.p_top_vs_maxnull); n_sig = int(sum(p < 0.05 for p in padj))
fig = plt.figure(figsize=(7.09, 3.54)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 180); ax.set_ylim(0, 90); ax.axis("off")
CA, CB, CP = ("#fdf3e1", "#c98a1f"), ("#e8f0fa", "#2f6db5"), ("#f0e9f6", "#7b4fa3"); X = [(4, 50), (61, 56), (124, 53)]; Y0, Hh = 6, 76
def box(x, w, col, title, k):
    ax.add_patch(FancyBboxPatch((x, Y0), w, Hh, boxstyle="round,pad=0,rounding_size=2.5", fc=col[0], ec=col[1], lw=0.9)); ax.add_patch(Circle((x + 5.5, Y0 + Hh - 6), 2.6, color=col[1])); ax.text(x + 5.5, Y0 + Hh - 6.1, str(k), color="white", fontsize=7.5, fontweight="bold", ha="center", va="center")
    ax.text(x + 10, Y0 + Hh - 6, title, color=col[1], fontsize=8, fontweight="bold", va="center")
def arrow(x0, x1, y=44, col="0.3"): ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>", mutation_scale=9, lw=1.1, color=col))
def check(x, y, col="#3a8f5c"): ax.plot([x - 1.6, x - 0.4, x + 1.9], [y, y - 1.5, y + 1.8], color=col, lw=1.6, solid_capstyle="round")
def cross(x, y, col="#b5483a"): ax.plot([x - 1.5, x + 1.5], [y - 1.5, y + 1.5], color=col, lw=1.6, solid_capstyle="round"); ax.plot([x - 1.5, x + 1.5], [y + 1.5, y - 1.5], color=col, lw=1.6, solid_capstyle="round")
def half(x, y, col="#c98a1f"): ax.add_patch(Circle((x, y), 1.7, fc="white", ec=col, lw=1.0)); ax.add_patch(matplotlib.patches.Wedge((x, y), 1.7, 90, 270, fc=col, ec=col, lw=0))
# 1 compendium
x, w = X[0]; box(x, w, CA, "Public compendium", 1); cx = x + 13
ax.add_patch(Rectangle((cx - 7, 50), 14, 14, fc="white", ec="none")); ax.add_patch(Ellipse((cx, 50), 14, 4, fc="white", ec=CA[1], lw=0.8)); ax.add_patch(Rectangle((cx - 6.95, 50), 13.9, 2.2, fc="white", ec="none"))
ax.plot([cx - 7, cx - 7], [50, 64], color=CA[1], lw=0.8); ax.plot([cx + 7, cx + 7], [50, 64], color=CA[1], lw=0.8); ax.add_patch(Ellipse((cx, 64), 14, 4, fc="white", ec=CA[1], lw=0.8))
for yy in (54.5, 59): ax.add_patch(matplotlib.patches.Arc((cx, yy), 14, 4, theta1=180, theta2=360, ec=CA[1], lw=0.6))
ax.text(x + 23, 61, "ENCODE", fontsize=7.5, fontweight="bold", va="center"); ax.text(x + 23, 56.5, f"{n_runs} knockdowns", fontsize=6.5, va="center"); ax.text(x + 23, 52.5, f"{n_rbp} RBPs, 2 cell lines", fontsize=6.5, va="center")
ax.text(x + 5, 42, "Two artefacts removed", fontsize=6.8, fontweight="bold", va="center")
for yy, t in ((36.5, "shared control"), (31, "control regression to the mean")): ax.text(x + 9, yy, t, fontsize=6.3, va="center", color="0.25"); cross(x + 6, yy)
bx = x + 8; ax.plot([bx, bx + 36], [12.5, 12.5], color="0.3", lw=0.6); sc = 7 / (rule_sig - 0.5)
for xi, v, lab, col in ((bx + 5, rule_sig, "before", "#c98a1f"), (bx + 21, rule_clean, "after", "#e5c98f")):
    ax.add_patch(Rectangle((xi, 12.5), 10, (v - 0.5) * sc + 0.4, fc=col, ec="none")); ax.text(xi + 5, 12.5 + (v - 0.5) * sc + 1.6, f"{v:.2f}", fontsize=6.3, ha="center", va="bottom"); ax.text(xi + 5, 10.2, lab, fontsize=5.8, ha="center", va="center", color="0.3")
ax.text(x + 5, 25.6, "accuracy of a zero-parameter rule", fontsize=5.8, color="0.3", va="center"); ax.text(bx + 37, 12.5, "0.50", fontsize=5.3, color="0.35", va="center")
# 2 forward
x, w = X[1]; box(x, w, CB, "Predict an unseen RBP", 2); ax.add_patch(Circle((x + 9, 60), 4.2, fc="white", ec=CB[1], lw=1.0)); ax.text(x + 9, 59.8, "?", fontsize=10, fontweight="bold", color=CB[1], ha="center", va="center")
ax.text(x + 16, 62.2, "never knocked down", fontsize=6.5, va="center"); ax.text(x + 16, 57.8, "in the training data", fontsize=6.5, va="center")
def ladder(x0, y, lo, hi, vals, refv, ctl, title, colr):
    L = 21; f = lambda v: x0 + (v - lo) / (hi - lo) * L; ax.plot([x0, x0 + L], [y, y], color="0.6", lw=0.7); ax.text(x0, y + 8.6, title, fontsize=6.6, fontweight="bold", va="center")
    ax.plot([f(ctl), f(ctl)], [y - 1.6, y + 1.6], color="0.35", lw=0.8); ax.text(f(ctl), y - 3.6, f"{ctl:.2f}", fontsize=5.3, ha="center", color="0.35"); ax.scatter([f(refv)], [y], marker="D", s=14, color="k", zorder=3, lw=0); ax.text(f(refv), y - 3.6, f"{refv:.2f}", fontsize=5.3, ha="center")
    for v, c in zip(vals, CELLS): ax.scatter([f(v)], [y], s=22, color=COL[c], zorder=4, edgecolor="white", lw=0.3)
    ax.text(x0, y + 4.4, " / ".join(f"{v:.2f}" for v in vals), fontsize=6.6, color=colr, va="center", fontweight="bold")
ladder(x + 5, 34, 0.48, 0.60, ev, ceil.mean(), 0.5, "Single events", "#b5483a"); ladder(x + 31, 34, 0.0, 0.70, r1, ref, perm, "Programmes", "#3a8f5c"); cross(x + 23.5, 43); check(x + 51.5, 43)
ax.text(x + 5, 25.2, "direction accuracy", fontsize=5.5, color="0.3", va="center"); ax.text(x + 31, 25.2, "programme-1 score, r", fontsize=5.5, color="0.3", va="center")
nx, ny = x + 8, 14; pts = np.array([[0, 0], [5, 3], [4, -3], [9, 0.5], [12, -2.5]]) + [nx, ny]
for a_, b_ in ((0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (3, 4)): ax.plot(pts[[a_, b_], 0], pts[[a_, b_], 1], color=CB[1], lw=0.6, zorder=1)
ax.scatter(pts[:, 0], pts[:, 1], s=16, color="white", edgecolor=CB[1], lw=0.8, zorder=2); ax.scatter(pts[3:4, 0], pts[3:4, 1], s=16, color=CB[1], zorder=3, lw=0)
ax.text(x + 24, 15.8, "what transfers:", fontsize=6.3, va="center", color="0.25"); ax.text(x + 24, 12, "interaction-network context", fontsize=6.3, va="center", fontweight="bold", color=CB[1])
# 3 inverse
x, w = X[2]; box(x, w, CP, "Name the perturbed RBP", 3)
for k, hgt in enumerate((5, -4, 6, -3, 2, -5)): ax.add_patch(Rectangle((x + 5 + k * 2.3, 60), 1.7, hgt * 0.8, fc=CP[1], ec="none"))
ax.plot([x + 4.5, x + 19.5], [60, 60], color="0.3", lw=0.5); ax.text(x + 22, 62.2, "splicing shift", fontsize=6.5, va="center"); ax.text(x + 22, 57.8, "in a tissue", fontsize=6.5, va="center"); ax.text(x + 5, 50.5, f"matched against {n_rbp} measured responses", fontsize=5.8, va="center", color="0.3")
rows = [(41.5, check, "ENCODE, across cell lines", f"top-5 hit rate {100 * top5:.0f} % vs {100 * exp5:.0f} % by chance"), (30, half, "Tumours, knockout cells", f"class level; rank fraction {(T.gene_rank / T.n_rbps).median():.2f} over {len(T)} tests"), (18.5, cross, "Premalignant tissue", f"{n_sig} of {len(padj)} adjusted analyses significant")]
for yy, sym, a_, b_ in rows: sym(x + 6.5, yy); ax.text(x + 11, yy + 1.9, a_, fontsize=6.6, fontweight="bold", va="center"); ax.text(x + 11, yy - 2.1, b_, fontsize=5.5, va="center", color="0.3")
arrow(55, 60); arrow(118, 123)
save(fig, "graphical_abstract")
