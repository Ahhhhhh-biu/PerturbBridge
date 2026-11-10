"""Correspondence of response programmes between cell lines: |cosine| between HepG2 and K562 programme loadings over the 21,077 shared events.
The matrix is parsed from the last block of 13_model/programme/programme_annotation.txt (written by 04_annotate_programmes.py)."""
from _common import *
L = open(f"{PB}/13_model/programme/programme_annotation.txt").read().split("== correspondence between cell lines")[1].strip().splitlines(); rows = [l.split() for l in L[2:12]]
M = np.array([[float(x) for x in r[1:]] for r in rows]); fig, ax = plt.subplots(figsize=(3.4, 3.1)); im = ax.imshow(M, cmap="Purples", vmin=0, vmax=0.8)
for i in range(10):
    for j in range(10): ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=4.8, color="w" if M[i, j] > 0.45 else "k")
ax.set_xticks(range(10)); ax.set_xticklabels([f"P{j + 1}" for j in range(10)], fontsize=6); ax.set_yticks(range(10)); ax.set_yticklabels([f"P{j + 1}" for j in range(10)], fontsize=6); ax.tick_params(length=0)
ax.set_xlabel("K562 programme"); ax.set_ylabel("HepG2 programme"); ax.set_title("Programme loadings, |cosine| between cell lines", loc="left", fontweight="bold", fontsize=7.5)
for sp in ax.spines.values(): sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.ax.tick_params(labelsize=6); fig.tight_layout(); save(fig, "fig_programme_correspondence")
