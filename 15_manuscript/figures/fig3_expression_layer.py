"""Figure 3: gene-expression output layer.
(a) measurement ceilings of the control-centred expression response: within-run split half (13_model/expr_within_run_<cell>.tsv), the same RBP in two
    independent experiments with the same or a different knockdown technology (expr_ceiling_cross_run_<cell>.tsv, one value per pair), and pairs of
    different RBPs (expr_zero_point_<cell>.tsv).
(b) held-out RBP-mean Pearson r of the full model, the shuffled-RBP control and the gene-only control (13_model/expr_vs_controls_by_rbp.tsv; mean and
    bootstrap 95 % CI over RBPs).
(c) named RBPs: full model against its shuffled-RBP control (same file)."""
import numpy as np, pandas as pd
from _style import *
rng = np.random.default_rng(0); fig = plt.figure(figsize=(7.09, 2.7)); gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 0.8, 1.35], wspace=0.5, left=0.085, right=0.99, top=0.88, bottom=0.27)
ax = fig.add_subplot(gs[0]); cats = ["within run,\nsplit half", "same RBP,\nsame\ntechnology", "same RBP,\ncross\ntechnology", "different\nRBPs"]
for i, c in enumerate(CELLS):
    w = pd.read_csv(f"{PB}/13_model/expr_within_run_{c}.tsv", sep="\t").r_split_half.dropna(); x = pd.read_csv(f"{PB}/13_model/expr_ceiling_cross_run_{c}.tsv", sep="\t"); x = x.drop_duplicates(["target", "r_all"]); same = x.assay_select == x.assay_test
    zp = pd.read_csv(f"{PB}/13_model/expr_zero_point_{c}.tsv", sep="\t").r_all.dropna(); groups = [w, x.r_all[same], x.r_all[~same], zp]
    for j, g in enumerate(groups):
        xx = j + (i - 0.5) * 0.36; ax.bar(xx, g.mean(), 0.34, color=COL[c], label=c if j == 0 else None)
        if len(g) <= 30: ax.scatter(xx + rng.uniform(-0.09, 0.09, len(g)), g, s=5, color="k", alpha=0.55, lw=0, zorder=3)
        ax.text(xx, -0.47, f"{len(g)}", ha="center", fontsize=5.2, color="0.35")
ax.text(-0.62, -0.47, "n", ha="center", fontsize=5.2, color="0.35", style="italic"); ax.axhline(0, color="k", lw=0.5); ax.set_ylim(-0.52, 0.62); ax.set_xticks(range(4)); ax.set_xticklabels(cats, fontsize=5.8); ax.set_ylabel("Pearson r of expression responses\nbetween two measurements", fontsize=6.5); ax.legend(fontsize=6, loc="upper right", handlelength=1.2); panel(ax, "a", "Measurement ceilings")
E = pd.read_csv(f"{PB}/13_model/expr_vs_controls_by_rbp.tsv", sep="\t"); ax = fig.add_subplot(gs[1]); mods = [("A_r_all", "full\nmodel"), ("shuf_r_all", "shuffled\nRBP"), ("none_r_all", "gene\nonly")]
for i, c in enumerate(CELLS):
    e = E[E.cell_line == c]
    for j, (m, _) in enumerate(mods):
        v = e[m].dropna().values; bs = v[rng.integers(0, len(v), (2000, len(v)))].mean(1); xx = j + (i - 0.5) * 0.36
        ax.bar(xx, v.mean(), 0.34, color=COL[c]); ax.errorbar(xx, v.mean(), yerr=[[v.mean() - np.percentile(bs, 2.5)], [np.percentile(bs, 97.5) - v.mean()]], color="k", lw=0.6, capsize=1.5, capthick=0.6)
ax.axhline(0, color="k", lw=0.5); ax.set_xticks(range(3)); ax.set_xticklabels([m[1] for m in mods], fontsize=5.8); ax.set_ylabel("Held-out RBP-mean Pearson r"); ax.set_ylim(-0.03, 0.065); panel(ax, "b", "Unseen RBPs")
ax = fig.add_subplot(gs[2]); names = ["SF3B1", "SF3A1", "U2AF2", "SRSF1", "PTBP1", "HNRNPC"]; k = 0; ticks = []; labs = []
for nme in names:
    for i, c in enumerate(CELLS):
        r = E[(E.cell_line == c) & (E.target == nme)]
        if len(r) == 0: continue
        r = r.iloc[0]; ax.plot([r.shuf_r_all, r.A_r_all], [k, k], color="0.6", lw=0.8, zorder=1); ax.scatter([r.shuf_r_all], [k], s=16, facecolor="white", edgecolor=COL[c], lw=0.9, zorder=2); ax.scatter([r.A_r_all], [k], s=16, color=COL[c], zorder=3, lw=0)
        ticks.append(k); labs.append(f"{nme}  {c}"); k += 1
    k += 0.5
ax.axvline(0, color="k", lw=0.5); ax.set_yticks(ticks); ax.set_yticklabels(labs, fontsize=5.6); ax.invert_yaxis(); ax.set_xlabel("Held-out Pearson r"); ax.set_xlim(-0.45, 0.85)
ax.scatter([], [], s=16, color="0.3", lw=0, label="full model"); ax.scatter([], [], s=16, facecolor="white", edgecolor="0.3", lw=0.9, label="shuffled RBP"); ax.legend(fontsize=5.8, loc="lower right", handletextpad=0.2); panel(ax, "c", "Named RBPs")
save(fig, "fig3_expression_layer")
