"""Figure 1: two artefacts that inflate apparent predictability, and the clean target.
(a) dPSI correlation between pairs of knockdown runs of different RBPs, sharing a control or not, before and after control-centring
    (means and SDs parsed from logs/control_centering.txt, written by 12_baselines/03_control_centering.py).
(b) direction accuracy on FDR-significant rows and on the clean target for the zero-parameter control-deviation rule and the v2 models
    (13_model/clean_target_summary.tsv, written by 13_model/07_clean_target.py; row-level, seed 0, held-out folds).
(c) schematic of the clean-target definition (thresholds as in 13_model/07_clean_target.py: each knockdown replicate deviates > 0.05 from the
    control-group consensus in the same direction and |control-centred dPSI| > 0.1)."""
import re, numpy as np, pandas as pd
from _style import *
from matplotlib.patches import FancyBboxPatch
fig = plt.figure(figsize=(7.09, 2.75)); gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.25, 1.2], wspace=0.42, left=0.085, right=0.99, top=0.88, bottom=0.25)
# ---- a: shared control
txt = open(f"{PB}/logs/control_centering.txt").read().split("\n"); vals = {}; cell = None; kind = None
for l in txt:
    m = re.match(r"(HepG2|K562) runs per control", l)
    if m: cell = m.group(1)
    for key in ("r_raw", "r_cc"):
        if l.startswith(key): kind = key
    t = l.replace("r_raw", "").replace("r_cc", "").split()
    if cell and kind and len(t) == 3 and t[0] in ("count", "mean", "std"): vals[(cell, kind, t[0])] = (float(t[1]), float(t[2]))     # (different control, same control)
ax = fig.add_subplot(gs[0]); groups = [("r_raw", 0, "different\ncontrol"), ("r_raw", 1, "same\ncontrol"), ("r_cc", 0, "different\ncontrol"), ("r_cc", 1, "same\ncontrol")]
for i, c in enumerate(CELLS):
    x = np.arange(4) + (i - 0.5) * 0.36; y = [vals[(c, k, "mean")][j] for k, j, _ in groups]; e = [vals[(c, k, "std")][j] for k, j, _ in groups]
    ax.bar(x, y, 0.34, color=COL[c], label=c, yerr=e, error_kw=dict(lw=0.6, capsize=1.5, capthick=0.6))
ax.axhline(0, color="k", lw=0.5); ax.set_xticks(range(4)); ax.set_xticklabels([g[2] for g in groups], fontsize=6); ax.set_ylabel("dPSI correlation between runs\nof different RBPs (mean ± SD)", fontsize=6.5); ax.set_ylim(-0.25, 0.68)
ax.text(0.5, 0.64, "raw dPSI", ha="center", fontsize=6.5); ax.text(2.5, 0.64, "control-centred", ha="center", fontsize=6.5); ax.axvline(1.5, color="0.7", lw=0.5, ls=":"); ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(0.52, 0.86), handlelength=1.2)
n = {c: [int(v) for v in vals[(c, "r_raw", "count")]] for c in CELLS}; ax.text(0.5, -0.31, f"run pairs, different / same control: {n['HepG2'][0]:,} / {n['HepG2'][1]} (HepG2), {n['K562'][0]:,} / {n['K562'][1]} (K562)", transform=ax.transAxes, fontsize=5.2, color="0.35", ha="center")
panel(ax, "a", "Shared-control artefact")
# ---- b: regression to the mean
S = pd.read_csv(f"{PB}/13_model/clean_target_summary.tsv", sep="\t").set_index(["cell_line", "set"]); ax = fig.add_subplot(gs[1])
preds = [("rule", "control-deviation rule\n(no parameters)", "#d08a1d"), ("event_only", "event-only model", "#9a9a9a"), ("v2_A", "full model (v2)", "#2f5d8a"), ("shuffled", "shuffled-RBP model", "#c9c9c9")]; sets = [("sig_FDR", "FDR-significant rows"), ("clean_0.05", "clean target")]
w = 0.1
for si, (s, sl) in enumerate(sets):
    for pi, (p, pl, pc) in enumerate(preds):
        for ci, c in enumerate(CELLS):
            x = si + (pi - 1.5) * 0.22 + (ci - 0.5) * w; ax.bar(x, S.loc[(c, s), p], w * 0.95, color=pc, hatch=None if ci == 0 else "////", edgecolor="white" if ci == 0 else "0.25", linewidth=0.3, label=pl if (si == 0 and ci == 0) else None)
ax.axhline(0.5, color="k", lw=0.5, ls="--"); ax.set_ylim(0.4, 0.9); ax.set_xticks([0, 1]); ax.set_xticklabels([s[1] for s in sets]); ax.set_ylabel("Direction accuracy, held-out rows")
h, l = ax.get_legend_handles_labels(); from matplotlib.patches import Patch
h += [Patch(facecolor="0.6", edgecolor="white"), Patch(facecolor="0.6", hatch="////", edgecolor="0.25", linewidth=0.3)]; l += ["HepG2 (plain)", "K562 (hatched)"]; ax.legend(h, l, fontsize=5.6, loc="upper right", bbox_to_anchor=(1.02, 1.0), handlelength=1.2, labelspacing=0.35)
ax.text(-0.33, S.loc[("HepG2", "sig_FDR"), "rule"] + 0.012, f"{S.loc[('HepG2', 'sig_FDR'), 'rule']:.2f} / {S.loc[('K562', 'sig_FDR'), 'rule']:.2f}", ha="center", fontsize=5.6, color="#d08a1d")
ax.text(1 - 0.33, S.loc[("HepG2", "clean_0.05"), "rule"] + 0.012, f"{S.loc[('HepG2', 'clean_0.05'), 'rule']:.2f} / {S.loc[('K562', 'clean_0.05'), 'rule']:.2f}", ha="center", fontsize=5.6, color="#d08a1d")
panel(ax, "b", "Control regression to the mean")
# ---- c: schematic
ax = fig.add_subplot(gs[2]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); panel(ax, "c", "Clean target")
x0 = 0.74; L = 0.50; rows = [(0.74, "clean", [0.87, 0.92], 0.75, "#3a8f5c", "both replicates shift\nthe same way"), (0.47, "not clean", [0.75, 0.72], 0.56, "#b5483a", "only the run's own\ncontrol is off"), (0.20, "not clean", [0.86, 0.60], 0.73, "#b5483a", "replicates disagree")]
ax.add_patch(plt.Rectangle((x0 - 0.05, 0.08), 0.10, 0.78, color="0.9", lw=0)); ax.plot([x0, x0], [0.08, 0.86], color="0.35", lw=0.8)
ax.text(1.0, 0.90, "control-group consensus ± 0.05", ha="right", fontsize=5.6, color="0.3")
for y, lab, reps, ctrl, colr, why in rows:
    ax.plot([L, 1.0], [y, y], color="0.75", lw=0.5); ax.scatter(reps, [y, y], s=22, color="#2f5d8a", zorder=3, lw=0); ax.scatter([ctrl], [y], s=22, marker="D", facecolor="white", edgecolor="0.2", lw=0.7, zorder=3)
    ax.text(0.0, y + 0.05, lab, fontsize=6.3, fontweight="bold", color=colr, ha="left", va="center"); ax.text(0.0, y - 0.045, why, fontsize=5.3, color="0.3", ha="left", va="center", linespacing=1.1)
ax.scatter([0.02], [0.0], s=22, color="#2f5d8a", lw=0, clip_on=False); ax.text(0.05, 0.0, "knockdown replicate", fontsize=5.6, va="center"); ax.scatter([0.02], [-0.07], s=22, marker="D", facecolor="white", edgecolor="0.2", lw=0.7, clip_on=False); ax.text(0.05, -0.07, "the run's own control", fontsize=5.6, va="center")
ax.annotate("", xy=(1.0, 0.0), xytext=(L, 0.0), arrowprops=dict(arrowstyle="-|>", lw=0.6, color="0.3"), annotation_clip=False); ax.text((L + 1) / 2, -0.07, "PSI of one event", ha="center", fontsize=6, color="0.3")
save(fig, "fig1_artefacts")
