"""All premalignant-cohort contrasts in one view. Left: observed top dictionary score (dot) against the median of the permutation-maximum null
(tick), with p_top. Right: number of RBPs at label-permutation p < 0.05 against the number expected by chance.
Sources: 14_inverse/h4_summary.tsv, h4new_summary.tsv and, if present, h4b2_summary.tsv (same columns as h4new_summary.tsv); null medians are
parsed from logs/h4*.txt ("permutation-max null median"). The 95th percentile of the null maximum is not stored in tables or logs and is not shown."""
from _common import *; import glob, os
def null_medians():
    out = {}
    for f in sorted(glob.glob(f"{PB}/logs/h4*.txt")):
        name = adj = None
        for l in open(f, errors="replace"):
            m = re.match(r"== ([\w.+-]+):", l)
            if m: name, adj = m.group(1), None
            m = re.match(r"-- (\w+):", l)
            if m: adj = m.group(1)
            m = re.search(r"permutation-max null median ([0-9.]+)", l)
            if m and name: out[f"{name}_{adj}" if adj else name] = float(m.group(1))
    return out
NM = null_medians(); rows = []
A = pd.read_csv(f"{PB}/14_inverse/h4_summary.tsv", sep="\t"); lab0 = {"GSE274628_fallopian": "GSE274628 fallopian tube: BRCA carrier vs non-carrier (41 vs 53)", "GSE166044_breast": "GSE166044 breast: later diagnosed vs healthy (15 vs 15)"}
for _, r in A.iterrows(): rows.append(dict(group=lab0.get(r.cohort, r.cohort), adj=r.adjustment, top=r.top_score, rbp=r.top_rbp, p=r.p_top_vs_maxnull, n05=r.n_rbp_perm_p05, exp=r.expected_p05, null=NM.get(f"{r.cohort}_{r.adjustment}"), bh=np.nan))
skipped = []
for f in ("h4new_summary.tsv", "h4b2_summary.tsv"):
    if not os.path.exists(f"{PB}/14_inverse/{f}"): continue
    for _, r in pd.read_csv(f"{PB}/14_inverse/{f}", sep="\t").iterrows():
        m = re.match(r"(GSE\d+)_(.*?)_(adjusted_plus_composition|adjusted|unadjusted|descriptive)$", r.contrast); gse, what, adj = m.groups() if m else ("", r.contrast, "")
        if not np.isfinite(r.top_score): skipped.append((gse, what, adj, int(r.n_query))); continue
        rows.append(dict(group=f"{gse} {what.replace('_', ' ')} ({int(r.n_group1)} vs {int(r.n_group0)})", adj=adj, top=r.top_score, rbp=r.top_rbp, p=r.p_top_vs_maxnull, n05=r.n_rbp_perm_p05, exp=r.expected_p05, null=NM.get(r.contrast), bh=r.get("n_bh_q10", np.nan)))
D = pd.DataFrame(rows); adjlab = {"unadjusted": "unadjusted", "adjusted": "adjusted", "adjusted_plus_composition": "adjusted + composition", "descriptive": "descriptive control"}; adjcol = {"unadjusted": "#bfc5cc", "adjusted": "#2f5d8a", "adjusted_plus_composition": "#123456", "descriptive": "#8e44ad"}
order = {"unadjusted": 0, "adjusted": 1, "adjusted_plus_composition": 2, "descriptive": 3}; groups = list(dict.fromkeys(D.group)); y = []; yy = 0; gpos = {}
for g in groups:
    gpos[g] = yy; yy += 1
    for i in D[D.group == g].sort_values("adj", key=lambda s: s.map(order)).index: y.append((i, yy)); yy += 1
    yy += 0.5
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.09, 0.62 + 0.145 * yy), gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.08}, sharey=True)
for g, gy in gpos.items(): a1.text(-0.02, gy, g, transform=a1.get_yaxis_transform(), ha="right", va="center", fontsize=6, fontweight="bold")
for i, yv in y:
    r = D.loc[i]; c = adjcol.get(r.adj, "k")
    if r.null is not None and np.isfinite(r.null): a1.plot([r.null, r.top], [yv, yv], color=c, lw=0.8, zorder=1); a1.plot([r.null], [yv], marker="|", ms=6, color="k", mew=0.9, zorder=2)
    a1.scatter([r.top], [yv], s=16, color=c, edgecolor="k", linewidth=0.3, zorder=3); a1.text(-0.02, yv, adjlab.get(r.adj, r.adj), transform=a1.get_yaxis_transform(), ha="right", va="center", fontsize=5.6)
    a1.text(max(r.top, r.null or 0) + 0.012, yv, f"{r.rbp}, p = {r.p:.2f}" if r.p >= 0.0095 else f"{r.rbp}, p = {r.p:.3f}", va="center", fontsize=5.4, fontweight="bold" if r.p < 0.05 else "normal")
    a2.barh(yv, r.n05, color=c, height=0.7); a2.text(max(r.n05, r.exp) + 1.5, yv, f"{int(r.n05)}" + (f"  (BH q < 0.1: {int(r.bh)})" if np.isfinite(r.bh) and r.bh > 0 else ""), va="center", fontsize=5.4)
a2.axvline(D.exp.median(), color="k", lw=0.6, ls="--"); a2.text(D.exp.median() + 1.5, -0.55, f"expected by chance ({D.exp.median():.0f})", ha="left", va="center", fontsize=5.6)
a1.set_ylim(yy - 0.3, -1.2); a1.set_yticks([]); a2.tick_params(left=False); a1.spines["left"].set_visible(False); a2.spines["left"].set_visible(False); a1.set_xlim(0.15, max(D.top.max(), 0.5) + 0.2); a2.set_xlim(0, D.n05.max() * 1.75)
a1.set_xlabel("Top dictionary score (Pearson r)\ndot: observed; tick: median of the permutation maximum"); a2.set_xlabel("Number of RBPs (of 478)")
a1.set_title("a  Best match vs chance (name, p_top)", loc="left", fontweight="bold", fontsize=7.5); a2.set_title("b  RBPs at permutation p < 0.05", loc="left", fontweight="bold", fontsize=7.5)
if skipped: fig.text(0.01, 0.005, "Not evaluable (fewer than 50 query events): " + "; ".join(sorted({f"{g} ({min(q for g2, w, a_, q in skipped if g2 == g)}–{max(q for g2, w, a_, q in skipped if g2 == g)} events, {sum(g2 == g for g2, *_ in skipped)} contrasts)" for g, *_ in skipped})), fontsize=5.6, va="bottom")
fig.subplots_adjust(left=0.375, right=0.985, top=1 - 0.28 / (0.5 + 0.145 * yy), bottom=0.62 / (0.5 + 0.145 * yy)); save(fig, "fig_precancer_forest"); print(D.to_string()); print("null medians found:", int(D.null.notna().sum()), "of", len(D), "| skipped", len(skipped))
