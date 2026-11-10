"""Second-stage H4 calls and cell composition (fallopian GSE274628, breast GSE166044).
a: RBP calls before and after adding all composition scores and age (14_inverse/h4rob_summary.tsv: n_prev, n_surviving; family-wise p before from
   h4ext_summary.tsv, after from h4rob_summary.tsv). b: t statistic of the fallopian splicing-layer calls under graded adjustment (h4adj_check.tsv).
c: composition-matched pairs without regression, paired t and sign-flip p (parsed from logs/h4_adjustment_check.txt, section 4; for CSNK1E the log
   gives only the maximum |t| over all RBPs, shown as |t|). d: sample-level activity against the ciliated-cell score for a lost call (RRP9) and for SRSF3
   (14_inverse/h4ext_activity_fallopian_splicing.parquet, 05_precancer_cohorts/GSE274628/cell_composition_scores.tsv; r computed here)."""
from _common import *; from matplotlib.gridspec import GridSpec
E = pd.read_csv(f"{PB}/14_inverse/h4ext_summary.tsv", sep="\t"); Rb = pd.read_csv(f"{PB}/14_inverse/h4rob_summary.tsv", sep="\t"); G = pd.read_csv(f"{PB}/14_inverse/h4adj_check.tsv", sep="\t")
fig = plt.figure(figsize=(7.09, 5.6)); gs = GridSpec(2, 4, figure=fig, hspace=0.62, wspace=0.75, left=0.075, right=0.985, top=0.94, bottom=0.085)
ax = fig.add_subplot(gs[0, :2]); lab = {"A_sample_level": "splicing", "D_expression": "expression"}; x = np.arange(len(Rb))
ax.bar(x - 0.24, Rb.n_prev, 0.4, color="#bfc5cc", label="batch, menopause, epithelial score"); ax.bar(x + 0.24, Rb.n_calls_full, 0.4, color="#2f5d8a", label="+ all composition scores and age")
for i, r in Rb.iterrows():
    pb = E[(E.cohort == r.cohort) & (E.analysis == r.layer)].p_familywise_maxt.iloc[0]; ax.text(i - 0.24, r.n_prev + 1.5, f"{int(r.n_prev)}\np = {pb:.3f}", ha="center", va="bottom", fontsize=5.0); ax.text(i + 0.24, r.n_calls_full + 1.5, f"{int(r.n_calls_full)}\np = {r.p_familywise:.3f}", ha="center", va="bottom", fontsize=5.0)
ax.set_xticks(x); ax.set_xticklabels([f"{r.cohort}\n{lab.get(r.layer, r.layer)}" for _, r in Rb.iterrows()], fontsize=6); ax.set_ylabel("RBP calls (BH q < 0.1)"); ax.set_ylim(0, Rb.n_prev.max() * 1.8); ax.legend(fontsize=5.6, loc="upper right", title="Covariates (p: family-wise, max |t|)", title_fontsize=5.8)
ax.text(0.98, 0.70, f"earlier calls surviving the full model: {int(Rb.n_surviving.sum())} of {int(Rb.n_prev.sum())}", transform=ax.transAxes, fontsize=5.6, va="top", ha="right")
ax.set_title("a  Calls before and after full composition adjustment", loc="left", fontweight="bold", fontsize=7.5)
ax = fig.add_subplot(gs[0, 2:]); cols = [c for c in G.columns if c.startswith("t_")]; pal = {"t_RRP9": "#c0392b", "t_EIF3D": "#e67e22", "t_EIF2S1": "#d4ac0d", "t_SRSF3": "#2f5d8a", "t_SUGP2": "#7f8c8d", "t_TAF15": "#8e44ad"}
for c in cols: ax.plot(range(len(G)), G[c], marker="o", ms=3, lw=1, color=pal.get(c, "k"), label=c[2:])
ax.axhline(0, color="k", lw=0.3, ls=":"); short = ["none", "batch +\nmenopause", "+ epithelial", "+ 4 cell-type\nscores", "all 10 scores\n+ age"]; ax.set_xticks(range(len(G))); ax.set_xticklabels([f"{s}\n(p = {p:.2f})" for s, p in zip(short[:len(G)], G.familywise_p)], fontsize=5.3)
ax.set_ylabel("t of carrier status on RBP activity"); ax.legend(fontsize=5.6, ncol=2, loc="center left", bbox_to_anchor=(0.02, 0.58), columnspacing=1, handlelength=1.5); ax.set_xlabel("Covariates in the model (family-wise p over 478 RBPs)", fontsize=6.2)
ax.set_title("b  Fallopian splicing-layer calls under graded adjustment", loc="left", fontweight="bold", fontsize=7.5)
L = open(f"{PB}/logs/h4_adjustment_check.txt").read().split("== 4.")[1]; mp = re.findall(r"(\w+)\s+mean paired diff ([+-][0-9.]+)\s+t ([+-][0-9.]+)\s+parametric p [0-9.]+\s+sign-flip p ([0-9.]+)", L); npairs = re.search(r"(\d+) carriers each matched", L).group(1); fw = re.search(r"max paired \|t\| ([0-9.]+), sign-flip p ([0-9.]+), top RBP (\w+)", L)
ax = fig.add_subplot(gs[1, :2]); names = [m[0] for m in mp]; tv = [float(m[2]) for m in mp]; pv = [float(m[3]) for m in mp]; ax.bar(range(len(mp)), tv, 0.62, color=["#2f5d8a" if p < 0.05 else "#bfc5cc" for p in pv])
for i, (t, p, m) in enumerate(zip(tv, pv, mp)): ax.text(i, t + (0.18 if t >= 0 else -0.18), f"diff {m[1]}\np = {p:.3f}", ha="center", va="bottom" if t >= 0 else "top", fontsize=5.2)
if fw: ax.bar(len(mp), float(fw.group(1)), 0.62, color="#2f5d8a", hatch="///", edgecolor="w", lw=0); ax.text(len(mp), float(fw.group(1)) + 0.18, f"largest |t| of 478\np = {float(fw.group(2)):.3f}", ha="center", va="bottom", fontsize=5.2); names.append(f"{fw.group(3)}\n(|t|)")
ax.axhline(0, color="k", lw=0.4); ax.set_xticks(range(len(names))); ax.set_xticklabels(names, fontsize=6); ax.set_ylabel(f"Paired t, carrier minus matched non-carrier\n({npairs} composition-matched pairs)"); ax.set_ylim(min(tv) - 1.6, (float(fw.group(1)) if fw else max(tv)) + 1.8)
ax.set_title("c  Matched pairs, no regression (sign-flip p)", loc="left", fontweight="bold", fontsize=7.5)
Aq = pd.read_parquet(f"{PB}/14_inverse/h4ext_activity_fallopian_splicing.parquet", columns=["RRP9", "SRSF3"]); C = pd.read_csv(f"{PB}/05_precancer_cohorts/GSE274628/cell_composition_scores.tsv", sep="\t").set_index("gsm").reindex(Aq.index); car = C.group.str.startswith("b1").values
for j, rbp in enumerate(("RRP9", "SRSF3")):
    ax = fig.add_subplot(gs[1, 2 + j]); xv = C.ft_ciliated.values; yv = Aq[rbp].values; ax.scatter(xv[~car], yv[~car], s=7, color="#bfc5cc", edgecolor="none", label=f"non-carrier (n = {(~car).sum()})"); ax.scatter(xv[car], yv[car], s=7, color="#c0392b", edgecolor="none", label=f"BRCA carrier (n = {car.sum()})")
    r_all = np.corrcoef(xv, yv)[0, 1]; r_nc = np.corrcoef(xv[~car], yv[~car])[0, 1]; ax.text(0.04, 0.04, f"r = {r_all:.2f} (all)\nr = {r_nc:.2f} (non-carriers)", transform=ax.transAxes, fontsize=5.6, va="bottom", bbox=dict(facecolor="w", edgecolor="none", alpha=0.8, pad=1))
    ax.set_xlabel("Ciliated-cell score (log2 CPM)", fontsize=6.2); ax.set_ylabel(f"{rbp} activity score", fontsize=6.2); ax.tick_params(labelsize=6); ax.set_title(("d  " if j == 0 else "") + f"{rbp}" + (" (call lost)" if j == 0 else " (kept by matching)"), loc="left", fontweight="bold" if j == 0 else "normal", fontsize=7.5 if j == 0 else 6.5)
    if j == 0: ax.legend(fontsize=5.2, loc="upper right", handletextpad=0.1, borderaxespad=0.1)
save(fig, "fig_composition_adjustment"); print(mp, fw.groups() if fw else None)
