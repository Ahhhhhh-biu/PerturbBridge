"""RBP x 10 programme score heatmap per cell line. Scores: 13_model/programme/programme_scores_<cell>.tsv (reference basis over all RBPs), z-scored per
programme and clipped at +/-4. Rows: RBPs with any |z| > 2.5 (the 60 with the largest max |z| if more), average-linkage clustering on Euclidean distance.
Side bar: complex class from Complex Portal (see _common.py)."""
from _common import *; from scipy.cluster.hierarchy import linkage, leaves_list
cls = rbp_classes(); fig, axes = plt.subplots(1, 2, figsize=(7.09, 8.0), gridspec_kw={"wspace": 0.42}); present = set()
for ax, c in zip(axes, ("HepG2", "K562")):
    S = pd.read_csv(f"{PB}/13_model/programme/programme_scores_{c}.tsv", sep="\t", index_col=0); Z = (S - S.mean()) / S.std(ddof=0); mx = Z.abs().max(1); sel = mx[mx > 2.5].sort_values(ascending=False).index[:60]
    Zs = Z.loc[sel]; Zs = Zs.iloc[leaves_list(linkage(Zs.clip(-4, 4).values, "average"))]; im = ax.imshow(Zs.clip(-4, 4).values, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(10)); ax.set_xticklabels([f"P{i + 1}" for i in range(10)], fontsize=6); ax.set_yticks(range(len(Zs))); ax.set_yticklabels(Zs.index, fontsize=4.8); ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    cl = [cls.get(t, "other") for t in Zs.index]; present |= set(cl)
    for i, k in enumerate(cl): ax.add_patch(plt.Rectangle((10 - 0.5 + 0.15, i - 0.5), 0.45, 1, color=CLASS_COL[k], clip_on=False, lw=0))
    ax.set_title(f"{'ab'[c == 'K562']}  {c}: {len(Zs)} of {int((mx > 2.5).sum())} RBPs with any |z| > 2.5", loc="left", fontweight="bold", fontsize=7.5)
fig.subplots_adjust(left=0.085, right=0.965, top=0.965, bottom=0.115)
cax = fig.add_axes([0.085, 0.052, 0.30, 0.011]); cb = fig.colorbar(im, cax=cax, orientation="horizontal"); cb.set_label("Programme score (z across RBPs, clipped at ±4)", fontsize=6.5); cb.ax.tick_params(labelsize=6)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=CLASS_COL[c], label=c) for c in CLASS_COL if c in present], loc="lower right", bbox_to_anchor=(0.965, 0.012), ncol=3, fontsize=5.8, title="Complex class (Complex Portal), side bar", title_fontsize=6.2, handlelength=1, columnspacing=1.2)
save(fig, "fig_programme_heatmap")
