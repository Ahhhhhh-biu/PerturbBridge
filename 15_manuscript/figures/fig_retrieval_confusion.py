"""Who is retrieved for whom: ENCODE cross-cell-line retrieval (14_inverse/h3_encode_retrieval.tsv, setting == "cross"; one row per query run with
the ordered top-5 dictionary RBPs). Queries are restricted to RBPs of high or mid reproducibility class (column repro_class of the query's own
cell line). Matrix: query RBP (rows) x dictionary RBP (columns) over the same RBP set, value = best position (1-5) of the dictionary RBP in the
query's top-5 list over the runs of that query RBP. Order: complex class (Complex Portal, _common.py), then name. Bars: what the top-1 hit is."""
from _common import *; from matplotlib.colors import ListedColormap, BoundaryNorm
R = pd.read_csv(f"{PB}/14_inverse/h3_encode_retrieval.tsv", sep="\t"); R = R[(R.setting == "cross") & R.repro_class.isin(["high", "mid"])].copy(); cls = rbp_classes(); same = shared_complex()
order = list(CLASS_COL); T = sorted(R.target.unique(), key=lambda t: (order.index(cls.get(t, "other")), t)); pos = {t: i for i, t in enumerate(T)}; M = np.full((len(T), len(T)), np.nan)
for t, top5 in zip(R.target, R.top5.fillna("")):
    for k, d in enumerate(top5.split(";")):
        if d in pos and not (M[pos[t], pos[d]] <= k + 1): M[pos[t], pos[d]] = k + 1
kind = np.where(R.top1 == R.target, "the query RBP itself", np.where([same(a, b) for a, b in zip(R.target, R.top1)], "another member of a shared complex", "other RBP"))
fig = plt.figure(figsize=(7.09, 7.55)); ax = fig.add_axes([0.105, 0.078, 0.66, 0.66 * 7.09 / 7.55]); cmap = ListedColormap(["#08306b", "#2171b5", "#4292c6", "#6baed6", "#9ecae1"]); im = ax.imshow(M, cmap=cmap, norm=BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5, 5.5], 5), interpolation="nearest")
fs = max(3.1, min(5, 320 / len(T))); ax.set_xticks(range(len(T))); ax.set_xticklabels(T, rotation=90, fontsize=fs); ax.set_yticks(range(len(T))); ax.set_yticklabels(T, fontsize=fs); ax.tick_params(length=0, pad=5); ax.tick_params(axis="x", pad=8)
ax.plot([-0.5, len(T) - 0.5], [-0.5, len(T) - 0.5], color="0.7", lw=0.3, zorder=0); present = set()
for i, t in enumerate(T):
    k = cls.get(t, "other"); present.add(k); ax.add_patch(plt.Rectangle((-2.0, i - 0.5), 1.2, 1, color=CLASS_COL[k], clip_on=False, lw=0)); ax.add_patch(plt.Rectangle((i - 0.5, len(T) - 0.5 + 0.8), 1, 1.2, color=CLASS_COL[k], clip_on=False, lw=0))
b = [i for i in range(1, len(T)) if cls.get(T[i], "other") != cls.get(T[i - 1], "other")]
for x in b: ax.axhline(x - 0.5, color="k", lw=0.3); ax.axvline(x - 0.5, color="k", lw=0.3)
ax.set_xlim(-0.5, len(T) - 0.5); ax.set_ylim(len(T) - 0.5, -0.5); ax.set_xlabel("Dictionary RBP (other cell line)", labelpad=4); ax.set_ylabel("Query RBP (knockdown profile)", labelpad=6); ax.set_title(f"b  Top-5 retrieval among {len(T)} RBPs of high or mid reproducibility ({len(R)} query runs)", loc="left", fontweight="bold", fontsize=7.5, pad=6)
for sp in ax.spines.values(): sp.set_visible(True); sp.set_linewidth(0.4)
cax = fig.add_axes([0.80, 0.53, 0.012, 0.12]); cb = fig.colorbar(im, cax=cax, ticks=[1, 2, 3, 4, 5]); cb.set_label("Best position in the top 5", fontsize=6); cb.ax.tick_params(labelsize=6, length=0); cb.ax.invert_yaxis()
ax2 = fig.add_axes([0.19, 0.825, 0.42, 0.095]); cats = ["the query RBP itself", "another member of a shared complex", "other RBP"]; colb = ["#08306b", "#c0392b", "#bfc5cc"]
for j, rc in enumerate(("high", "mid")):
    v = [np.mean(kind[(R.repro_class == rc).values] == c_) for c_ in cats]; left = 0
    for x, cc in zip(v, colb): ax2.barh(j, x, left=left, color=cc, height=0.62); ax2.text(left + x / 2, j, f"{x * 100:.0f}", ha="center", va="center", fontsize=5.5, color="w" if cc != "#bfc5cc" else "k") if x > 0.06 else None; left += x
ax2.set_yticks([0, 1]); ax2.set_yticklabels([f"high (n = {(R.repro_class == 'high').sum()})", f"mid (n = {(R.repro_class == 'mid').sum()})"], fontsize=6); ax2.set_xlim(0, 1); ax2.set_xticks([0, 0.5, 1]); ax2.set_xticklabels(["0", "50", "100"]); ax2.set_xlabel("Top-1 hit is … (% of query runs)", fontsize=6.5); ax2.invert_yaxis()
ax2.set_title("a  Identity of the top-ranked dictionary RBP", loc="left", fontweight="bold", fontsize=7.5); from matplotlib.patches import Patch
ax2.legend(handles=[Patch(color=cc, label=c_) for c_, cc in zip(cats, colb)], loc="center left", bbox_to_anchor=(1.04, 0.5), fontsize=5.8, handlelength=1)
fig.legend(handles=[Patch(color=CLASS_COL[c_], label=c_) for c_ in CLASS_COL if c_ in present], loc="upper left", bbox_to_anchor=(0.79, 0.46), fontsize=5.6, title="Complex class (side bars)", title_fontsize=6, handlelength=1)
save(fig, "fig_retrieval_confusion"); print("queries", len(R), "RBPs", len(T), pd.Series(kind).value_counts().to_dict())
