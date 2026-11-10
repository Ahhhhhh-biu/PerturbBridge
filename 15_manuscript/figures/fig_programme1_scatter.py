"""Predicted versus measured programme-1 score of held-out RBPs for the selected model (sequence/text + network, ridge).
Source: 13_model/programme/per_rbp_final.tsv (set seq_net, method ridge; p1 predicted, t1 measured, pooled over the five folds); correlations printed
are read from results_summary_final.tsv (r1, rho1, r1_trim). Both axes are multiplied by -1 if needed so that SF3B1 has a positive measured score."""
from _common import *
R = pd.read_csv(f"{PB}/13_model/programme/per_rbp_final.tsv", sep="\t"); R = R[(R.set == "seq_net") & (R.method == "ridge")]; S = pd.read_csv(f"{PB}/13_model/programme/results_summary_final.tsv", sep="\t")
fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.7)); tops = []
for ax, c in zip(axes, ("HepG2", "K562")):
    g = R[R.cell == c].copy(); sg = np.sign(g.t1[g.target == "SF3B1"].iloc[0]) if (g.target == "SF3B1").any() else 1; g["p"] = g.p1 * sg; g["t"] = g.t1 * sg
    for k in ("low", "mid", "high"): h = g[g.repro == k]; ax.scatter(h.t, h.p, s=7 if k == "low" else 13, color=REPRO_COL[k], edgecolor="none" if k == "low" else "k", linewidth=0.3, label=f"{k} reproducibility (n = {len(h)})", zorder={"low": 1, "mid": 2, "high": 3}[k])
    top = g.reindex(g.t.abs().sort_values(ascending=False).index).head(12); ax.set_ylim(g.p.min() - 0.08 * np.ptp(g.p), g.p.max() + 0.42 * np.ptp(g.p))
    s = S[(S.cell == c) & (S.set == "seq_net") & (S.method == "ridge")].iloc[0]
    ax.text(0.03, 0.97, f"Pearson r = {s.r1:.2f}\nSpearman ρ = {s.rho1:.2f}\nPearson without top 5 % = {s.r1_trim:.2f}\npermuted descriptors (95th pct): r = {s.r1_perm_p95:.2f}", transform=ax.transAxes, va="top", fontsize=5.8)
    ax.axhline(0, color="k", lw=0.3, ls=":"); ax.axvline(0, color="k", lw=0.3, ls=":"); ax.set_xlabel("Measured programme-1 score (held-out RBP)"); ax.set_ylabel("Predicted programme-1 score")
    ax.set_title(f"{'ab'[c == 'K562']}  {c} (n = {len(g)} held-out RBPs)", loc="left", fontweight="bold", fontsize=7.5); ax.legend(fontsize=5.6, loc="upper right", handletextpad=0.2); tops.append((ax, top))
fig.tight_layout()
for ax, top in tops: repel_labels(ax, top.t.values, top.p.values, top.target.values)
save(fig, "fig_programme1_scatter")
