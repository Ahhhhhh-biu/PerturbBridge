"""Ranked scores of all dictionary RBPs for the HeyA8 LARP1 knockout splicing shift (14_inverse/heya8_ranking.tsv: score = Pearson r with the
dictionary profile, perm_p = label-permutation p). Highlighted: the top hits named in the manuscript and the pre-registered comparators LARP4 and LARP7.
In-house validation material: outputs carry the prefix inhouse_ (not for the public repository)."""
from _common import *
D = pd.read_csv(f"{PB}/14_inverse/heya8_ranking.tsv", sep="\t").sort_values("rank"); hits = ["PCBP2", "SF3B1", "PRPF8", "SF3A1", "SF3B4"]; comp = ["LARP4", "LARP7"]
fig, ax = plt.subplots(figsize=(4.6, 3.3)); ax.scatter(D["rank"], D.score, s=4, color="#bfc5cc", edgecolor="none", zorder=1, label=f"dictionary RBPs (n = {len(D)})")
h = D[D.rbp.isin(hits)]; c = D[D.rbp.isin(comp)]; ax.scatter(h["rank"], h.score, s=18, color="#c0392b", edgecolor="k", linewidth=0.3, zorder=3, label="core spliceosome / top hits"); ax.scatter(c["rank"], c.score, s=18, color="#2471a3", edgecolor="k", linewidth=0.3, zorder=3, label="pre-registered comparators (LARP family)")
ax.axhline(0, color="k", lw=0.3, ls=":"); ax.set_xscale("log"); ax.set_xlim(0.8, len(D) * 1.15); ax.set_xlabel(f"Rank among {len(D)} dictionary RBPs (log scale)"); ax.set_ylabel("Similarity to the LARP1-knockout splicing shift (Pearson r)")
ax.set_title("HeyA8 LARP1 knockout: which knockdown does it resemble?", loc="left", fontweight="bold", fontsize=7.5); ax.legend(fontsize=5.6, loc="lower left"); fig.tight_layout()
lab = lambda r: f"{r.rbp}  (rank {int(r['rank'])}, p = {r.perm_p:.3f})" if "perm_p" in D and np.isfinite(r.perm_p) else f"{r.rbp}  (rank {int(r['rank'])})"
repel_labels(ax, h["rank"].values, h.score.values, [lab(r) for _, r in h.iterrows()], fontsize=5.4, dx=5)
for (_, r), off in zip(c.sort_values("rank").iterrows(), ((-70, -26), (-22, -40))): ax.annotate(lab(r), (r["rank"], r.score), xytext=off, textcoords="offset points", fontsize=5.4, arrowprops=dict(arrowstyle="-", lw=0.3, color="0.4", shrinkB=2))
sub = pd.concat([h, c])
fig.savefig(f"{FIG}/inhouse_fig_heya8_ranking.png"); fig.savefig(f"{FIG}/inhouse_fig_heya8_ranking.pdf"); print("saved", f"{FIG}/inhouse_fig_heya8_ranking.png"); print(sub[["rbp", "rank", "score", "perm_p"]].to_string(index=False))
