"""Figure 4 draft: inverse mapping. A: cross-cell retrieval rank of the true RBP by reproducibility class (ENCODE); B: same- vs cross-technology
dictionaries (reagent control); C: TCGA mutant-vs-wild-type tests (gene rank, if available); D: precancer cohorts top RBPs with permutation null
(if available). Output 14_inverse/fig_inverse_mapping.png/.pdf (300 dpi, Arial)."""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib import font_manager
PB = "/well/blagden/users/eic789/PerturbBridge"; O = f"{PB}/14_inverse"
for f in ("Arial.ttf", "Arial-Bold.ttf"): font_manager.fontManager.addfont(f"/well/blagden/users/eic789/.fonts/{f}")
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6, "savefig.dpi": 300, "pdf.fonttype": 42})
D = pd.read_csv(f"{O}/h3_encode_retrieval.tsv", sep="\t"); X = pd.read_csv(f"{O}/h3_cross_technology.tsv", sep="\t")
fig, ax = plt.subplots(2, 2, figsize=(7.2, 6.0))
a = ax[0, 0]; cr = D[D.setting == "cross"]; order = ["high", "mid", "low"]; cols = {"high": "#1f77b4", "mid": "#ff7f0e", "low": "#7f7f7f"}
for i, c in enumerate(order):
    v = cr[cr.repro_class == c].rank_frac.values; a.scatter(np.full(len(v), i) + rng.uniform(-0.25, 0.25, len(v)) if (rng := np.random.default_rng(0)) is not None else 0, v, s=5, c=cols[c], lw=0, alpha=0.7)
    a.plot([i - 0.3, i + 0.3], [np.median(v)] * 2, color="k", lw=1)
a.axhline(0.5, color="k", lw=0.5, ls=":"); a.set_yscale("log"); a.set_ylim(0.001, 1.05); a.set_xticks(range(3)); a.set_xticklabels([f"{c}\nn={int((cr.repro_class == c).sum())}" for c in order]); a.set_xlabel("RBP reproducibility class (query cell line)"); a.set_ylabel("rank of true RBP / number of candidates")
a.set_title("Cross-cell-line retrieval (ENCODE)", fontsize=8, loc="left"); a.text(2.4, 0.5, "random", fontsize=6, va="bottom")
a = ax[0, 1]; S = X.groupby(["query_assay", "dict_assay"]).agg(top5=("rank", lambda x: (x <= 5).mean()), top10=("rank", lambda x: (x <= 10).mean()), n=("run", "size")).reset_index(); S["pair"] = S.query_assay + "→" + S.dict_assay
x = np.arange(len(S)); a.bar(x - 0.2, S.top5, 0.4, color="#1f77b4", label="top-5"); a.bar(x + 0.2, S.top10, 0.4, color="#aec7e8", label="top-10"); a.axhline(0.022, color="#d62728", lw=0.6, ls="--", label="random top-10")
a.set_xticks(x); a.set_xticklabels([f"{p}\nn={n}" for p, n in zip(S.pair, S.n)], fontsize=5.5, rotation=15); a.set_ylabel("hit rate"); a.set_title("Same- vs cross-technology dictionary", fontsize=8, loc="left"); a.legend(frameon=False, fontsize=6)
a = ax[1, 0]
if os.path.exists(f"{O}/h3_tcga_tests.tsv"):
    T = pd.read_csv(f"{O}/h3_tcga_tests.tsv", sep="\t"); T = T[T.gene_rank.notna()] if "gene_rank" in T else T.iloc[0:0]
    if len(T):
        T["frac"] = T.gene_rank / T.n_rbps; T["cls"] = np.where(T.test.str.contains("hotspot"), "splicing-factor hotspot", "truncating"); rng2 = np.random.default_rng(1)
        for j, (c, col) in enumerate((("splicing-factor hotspot", "#d62728"), ("truncating", "#7f7f7f"))):
            v = T[T.cls == c]; a.scatter(np.full(len(v), j) + rng2.uniform(-0.2, 0.2, len(v)), v.frac, s=12, c=col, lw=0, alpha=0.8)
            for _, r in v[v.frac <= 0.02].iterrows(): a.annotate(r.test.replace("TCGA-", "").replace("_hotspot", "").replace("_truncating", "") + f" ({int(r.n_case)})", (j, r.frac), xytext=(6, 0), textcoords="offset points", fontsize=6, va="center")
        a.axhline(0.5, color="k", lw=0.5, ls=":"); a.set_yscale("log"); a.set_ylim(0.0015, 1.1); a.set_xlim(-0.6, 1.9); a.set_xticks([0, 1]); a.set_xticklabels([f"SF3B1 / U2AF1 hotspot\nn={int((T.cls != 'truncating').sum())}", f"truncating, any RBP\nn={int((T.cls == 'truncating').sum())}"])
        a.set_ylabel("rank of mutated RBP / number of candidates"); a.text(1.55, 0.5, "random", fontsize=6, va="bottom")
    a.set_title("TCGA: mutant vs wild-type tumours (median 0.44)", fontsize=8, loc="left")
else: a.text(0.5, 0.5, "TCGA pending", ha="center", transform=a.transAxes); a.set_axis_off()
a = ax[1, 1]; k = 0
for coh, lab in (("GSE274628_fallopian", "fallopian tube, BRCA carrier vs non-carrier"), ("GSE166044_breast", "normal breast, susceptible vs healthy")):
    f = f"{O}/h4_{coh}_adjusted_ranking.tsv"
    if not os.path.exists(f): continue
    r = pd.read_csv(f, sep="\t").head(8); y = np.arange(len(r)) + k; a.barh(y, r.z_vs_null, color="#bbbbbb")
    for yy, nm in zip(y, r.rbp): a.text(-0.1, yy, nm, ha="right", va="center", fontsize=6)
    sm = pd.read_csv(f"{O}/h4_summary.tsv", sep="\t"); pt = sm[(sm.cohort == coh) & (sm.adjustment == "adjusted")].p_top_vs_maxnull.iloc[0]
    a.text(0.02, k + len(r) - 0.5, f"{lab}: top score vs permutation-max null p = {pt:.2f}", fontsize=6, style="italic"); k += len(r) + 1.5
if k: a.set_yticks([]); a.set_xlabel("z vs group-permutation null"); a.axvline(0, color="k", lw=0.5); a.set_title("Precancer cohorts: top RBPs (none beyond null)", fontsize=8, loc="left")
else: a.text(0.5, 0.5, "cohorts pending", ha="center", transform=a.transAxes); a.set_axis_off()
for i, (r_, c_) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]): ax[r_, c_].text(-0.12, 1.06, "ABCD"[i], transform=ax[r_, c_].transAxes, fontsize=10, fontweight="bold")
plt.tight_layout(); plt.savefig(f"{O}/fig_inverse_mapping.png"); plt.savefig(f"{O}/fig_inverse_mapping.pdf"); print("fig", f"{O}/fig_inverse_mapping.png")
