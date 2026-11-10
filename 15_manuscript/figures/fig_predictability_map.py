"""Predictability map: one point per RBP and cell line. x = the RBP's own reproducibility (within-run direction agreement, agree_0_1), y = held-out
direction accuracy on clean rows of the v3 deep model; bin means of v3, of its shuffled-RBP control and of the event-only control.
Source: 13_model/benchmark/master_by_rbp.tsv (written by 13_model/benchmark/03_benchmark_metrics.py)."""
import numpy as np, pandas as pd
from _style import *
M = pd.read_csv(f"{PB}/13_model/benchmark/master_by_rbp.tsv", sep="\t"); fig, axs = plt.subplots(1, 2, figsize=(7.09, 3.3), sharey=True); fig.subplots_adjust(left=0.075, right=0.99, top=0.9, bottom=0.14, wspace=0.06)
LAB = {"HepG2": ["U2AF2", "U2AF1", "SF3A3", "SF3B1", "SF3A1", "WBP11", "PUF60", "SRSF1", "PTBP1", "HNRNPC", "RBM22", "HNRNPK", "PRPF8"], "K562": ["SF3B3", "SF3B1", "SF3A1", "U2AF2", "U2AF1", "AQR", "CDC5L", "EFTUD2", "SRSF1", "PTBP1", "HNRNPC", "PRPF8", "SF3B4", "WBP11", "PUF60"]}
BINS = [0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 1.00]
def running(d, col="sign_acc"):
    mid = []; val = []
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        q = d[(d.agree_0_1 >= lo) & (d.agree_0_1 < hi)]
        if len(q) >= 5: mid.append(q.agree_0_1.mean()); val.append(q[col].mean())
    return mid, val
def place(ax, pts):
    """Greedy label placement: try offsets around each point, keep the first whose box hits no earlier box and no labelled point."""
    boxes = []; H = 0.026
    for nme, x, y in sorted(pts, key=lambda t: -t[2]):
        W = 0.0085 * len(nme) + 0.006; cands = [(0.014, 0), (0.014, 0.03), (0.014, -0.03), (-0.014 - W, 0), (-0.014 - W, 0.03), (-0.014 - W, -0.03), (0.02, 0.06), (0.02, -0.06), (-0.02 - W, 0.06), (-0.02 - W, -0.06), (0.03, 0.09), (-0.03 - W, 0.09)]
        for dx, dy in cands:
            bx = (x + dx, y + dy - H / 2, x + dx + W, y + dy + H / 2)
            if bx[2] > 0.95 or bx[0] < 0.42 or bx[3] > 0.995: continue
            if all(bx[2] < o[0] or bx[0] > o[2] or bx[3] < o[1] or bx[1] > o[3] for o in boxes) and all(not (bx[0] - 0.004 < px < bx[2] + 0.004 and bx[1] - 0.008 < py < bx[3] + 0.008) for _, px, py in pts if (px, py) != (x, y)): break
        boxes.append(bx); ax.annotate(nme, (x, y), (bx[0], y + dy), fontsize=5.4, va="center", ha="left", zorder=6, bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.8), arrowprops=dict(arrowstyle="-", lw=0.3, color="0.3", shrinkA=0, shrinkB=1.5, relpos=(0 if dx > 0 else 1, 0.5)))
for ax, (i, c) in zip(axs, enumerate(CELLS)):
    v = M[(M.cell_line == c) & (M.model == "v3 deep")].dropna(subset=["agree_0_1", "sign_acc"]); sh = M[(M.cell_line == c) & (M.model == "v3 shuffled RBP")].dropna(subset=["agree_0_1", "sign_acc"]); eo = M[(M.cell_line == c) & (M.model == "v3 event-only")].dropna(subset=["agree_0_1", "sign_acc"])
    ax.scatter(v.agree_0_1, v.sign_acc, s=7, color=COL[c], alpha=0.75, lw=0, zorder=2); ax.axhline(0.5, color="k", lw=0.5, ls="--"); ax.axvline(0.55, color="0.75", lw=0.5, ls=":"); ax.axvline(0.7, color="0.75", lw=0.5, ls=":")
    for d, colr, lab, ls in ((v, "#1b2f45", "v3 deep", "-"), (sh, "#c0392b", "shuffled-RBP control", "-"), (eo, "#7f7f7f", "event-only control", "--")): ax.plot(*running(d), color=colr, lw=1.2, ls=ls, label=lab, zorder=3)
    ax.plot([0.42, 0.95], [0.42, 0.95], color="0.6", lw=0.5, label="accuracy = reproducibility")
    pts = [(nme, v[v.target == nme].agree_0_1.iloc[0], v[v.target == nme].sign_acc.iloc[0]) for nme in LAB[c] if (v.target == nme).any()]; ax.scatter([q[1] for q in pts], [q[2] for q in pts], s=9, facecolor="none", edgecolor="k", lw=0.5, zorder=4); place(ax, pts)
    ax.set_xlim(0.42, 0.95); ax.set_ylim(0.2, 1.0); ax.set_title(f"{'ab'[i]}  {c} ({len(v)} RBPs)", loc="left", fontweight="bold", fontsize=7.5)
    ax.text(0.485, 0.215, "not reproducible", fontsize=5.5, color="0.45", ha="center"); ax.text(0.83, 0.215, "reproducible", fontsize=5.5, color="0.45", ha="center")
axs[0].set_ylabel("Held-out direction accuracy on clean events"); axs[0].legend(fontsize=6, loc="upper left", handlelength=1.6); fig.supxlabel("Reproducibility of the RBP's own knockdown (within-run direction agreement); lines are bin means", fontsize=7, y=0.02)
save(fig, "fig_predictability_map")
