"""Two-dimensional map of RBPs in programme space: PCA (first two components) of the z-scored 10 programme scores per cell line
(13_model/programme/programme_scores_<cell>.tsv), coloured by complex class (Complex Portal, see _common.py)."""
from _common import *
cls = rbp_classes(); import matplotlib.patheffects as pe
fig, axes = plt.subplots(1, 2, figsize=(7.09, 4.3)); present = set(); labs = []
for ax, c in zip(axes, ("HepG2", "K562")):
    S = pd.read_csv(f"{PB}/13_model/programme/programme_scores_{c}.tsv", sep="\t", index_col=0); Z = ((S - S.mean()) / S.std(ddof=0)).values; U, s, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False); XY = U[:, :2] * s[:2]; ve = s ** 2 / (s ** 2).sum()
    if XY[list(S.index).index("SF3B1"), 0] < 0: XY[:, 0] *= -1
    k = np.array([cls.get(t, "other") for t in S.index]); present |= set(k); o = k == "other"; ax.scatter(XY[o, 0], XY[o, 1], s=5, color=CLASS_COL["other"], edgecolor="none", zorder=1)
    for kk in CLASS_COL:
        if kk != "other" and (k == kk).any(): ax.scatter(XY[k == kk, 0], XY[k == kk, 1], s=16, color=CLASS_COL[kk], edgecolor="k", linewidth=0.3, zorder=3)
    d = np.sqrt((XY ** 2).sum(1)); lab = set(np.array(S.index)[np.argsort(-d)[:18]]) | {t for t in ("SF3B1", "U2AF1", "U2AF2") if t in S.index}
    labs.append((ax, [list(S.index).index(t) for t in lab], XY, list(S.index)))
    ax.set_xlabel(f"PC1 of programme scores ({ve[0] * 100:.0f} %)"); ax.set_ylabel(f"PC2 ({ve[1] * 100:.0f} %)"); ax.set_title(f"{'ab'[c == 'K562']}  {c} ({len(S)} RBPs)", loc="left", fontweight="bold", fontsize=7.5)
fig.tight_layout(rect=(0, 0.11, 1, 1))
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=CLASS_COL[c], label=c) for c in CLASS_COL if c in present], loc="lower center", ncol=5, fontsize=5.8, title="Complex class (Complex Portal); labels: the 18 RBPs farthest from the origin", title_fontsize=6.2, handlelength=1, columnspacing=1.2)
for ax, ii, XY, names in labs: repel_labels(ax, XY[ii, 0], XY[ii, 1], [names[i] for i in ii], pad=1.35, path_effects=[pe.withStroke(linewidth=1.3, foreground="w")])
save(fig, "fig_programme_map")
