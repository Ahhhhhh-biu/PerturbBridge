"""Formal reproducibility-stratified report of the v3 forward model (clean target, no control-PSI input) against its controls and
measurement ceilings. Inputs: 13_model/v3_vs_controls_by_rbp.tsv (per RBP: v3 A/B/C, v3 event-only, v3 shuffled-RBP, v2 A, within-run
replicate direction agreement agree_0_1), clean_ceiling_cross_run.tsv (same RBP, independent run), clean_zero_point.tsv (different RBPs).
Outputs (13_model/): stratified_results.tsv  = per cell line x reproducibility bin x model: n RBPs, mean sign accuracy, bootstrap 95% CI,
                                               paired difference to shuffled-RBP control with CI
                     stratified_named_rbps.tsv = core spliceosome and sequence-specific RBPs, model vs controls vs own reproducibility
                     fig_stratified_v3.png (300 dpi, Arial) = A: per-RBP scatter reproducibility vs v3 A; B: bin means with CI for all models,
                                               with clean-set zero point and same-RBP ceiling; C: named RBPs"""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib import font_manager
PB = "/well/blagden/users/eic789/PerturbBridge"; M = f"{PB}/13_model"; rng = np.random.default_rng(0)
for f in ("Arial.ttf", "Arial-Bold.ttf", "Arial-Italic.ttf"): font_manager.fontManager.addfont(f"/well/blagden/users/eic789/.fonts/{f}")
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "savefig.dpi": 300, "pdf.fonttype": 42})
d = pd.read_csv(f"{M}/v3_vs_controls_by_rbp.tsv", sep="\t"); ceil = pd.read_csv(f"{M}/clean_ceiling_cross_run.tsv", sep="\t"); zero = pd.read_csv(f"{M}/clean_zero_point.tsv", sep="\t")
models = {"v3_A": "v3 full (tier A)", "v3_B": "v3 no eCLIP (tier B)", "v3_C": "v3 no eCLIP, no motif (tier C)", "v3_none": "v3 event-only", "v3_shuf": "v3 shuffled RBP", "v2_A": "v2 full (old target)"}
bins = [(0, 0.55, "<0.55"), (0.55, 0.6, "0.55-0.6"), (0.6, 0.7, "0.6-0.7"), (0.7, 1.01, ">0.7")]
def boot_mean(x, n=2000):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 2: return np.nan, np.nan, np.nan
    b = rng.choice(x, (n, len(x))).mean(1); return x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5)
out = []
for cl, g in d.groupby("cell_line"):
    for lo, hi, lab in bins + [(0, 1.01, "all"), (0.6, 1.01, ">0.6")]:
        s = g[(g.agree_0_1 >= lo) & (g.agree_0_1 < hi)]
        for m, name in models.items():
            mu, l, h = boot_mean(s[m]); dm, dl, dh = boot_mean(s[m] - s["v3_shuf"])
            out.append(dict(cell_line=cl, repro_bin=lab, model=m, model_label=name, n_rbp=len(s), sign_acc=mu, ci_lo=l, ci_hi=h, diff_vs_shuffled=dm, diff_ci_lo=dl, diff_ci_hi=dh))
S = pd.DataFrame(out); S.to_csv(f"{M}/stratified_results.tsv", sep="\t", index=False, float_format="%.4f")
named = ["SF3B1", "SF3A1", "SF3B3", "U2AF1", "U2AF2", "PRPF8", "AQR", "PTBP1", "SRSF1", "HNRNPC", "QKI", "RBFOX2", "LARP4"]
N = d[d.target.isin(named)][["cell_line", "target", "agree_0_1", "v3_A", "v3_B", "v3_C", "v3_none", "v3_shuf", "v2_A", "n_sig"]].copy()
N["group"] = np.where(N.target.isin(["SF3B1", "SF3A1", "SF3B3", "U2AF1", "U2AF2", "PRPF8", "AQR"]), "core spliceosome", "sequence-specific")
N = N.sort_values(["cell_line", "group", "target"]); N.to_csv(f"{M}/stratified_named_rbps.tsv", sep="\t", index=False, float_format="%.3f")
ref = {cl: dict(zero=zero[zero.cell_line == cl].sign_agree.mean(), ceil=ceil[ceil.cell_line == cl].sign_agree.mean(), ceil_both=ceil[ceil.cell_line == cl].sign_agree_both_clean.mean()) for cl in ("HepG2", "K562")}
# figure
fig, ax = plt.subplots(3, 2, figsize=(7.2, 8.2))
for j, cl in enumerate(("HepG2", "K562")):
    g = d[d.cell_line == cl]; a = ax[0, j]
    a.scatter(g.agree_0_1, g.v3_shuf, s=6, c="#bbbbbb", lw=0, label="shuffled-RBP control"); a.scatter(g.agree_0_1, g.v3_A, s=7, c="#1f77b4", lw=0, label="v3 full model")
    for x in (0.55, 0.6, 0.7): a.axvline(x, color="grey", lw=0.4, ls=":")
    a.axhline(0.5, color="k", lw=0.4); a.set_xlim(0.38, 0.95); a.set_ylim(0.1, 1.0); a.set_xlabel("RBP's own replicate direction agreement (within run)"); a.set_ylabel("held-out direction accuracy, clean set")
    a.set_title(f"{cl}: {len(g)} RBPs", fontsize=8, loc="left"); a.legend(frameon=False, fontsize=6, loc="upper left")
    for tg in ("SF3B1", "SF3A1", "U2AF2", "PTBP1", "SRSF1", "HNRNPC"):
        r = g[g.target == tg]
        if len(r): a.annotate(tg, (r.agree_0_1.iloc[0], r.v3_A.iloc[0]), fontsize=5.5, xytext=(2, 2), textcoords="offset points")
    a = ax[1, j]; cols = {"v3_A": "#1f77b4", "v3_none": "#2ca02c", "v3_shuf": "#7f7f7f", "v2_A": "#ff7f0e"}; w = 0.2
    for k, (m, c) in enumerate(cols.items()):
        s = S[(S.cell_line == cl) & (S.model == m) & S.repro_bin.isin([b[2] for b in bins])].set_index("repro_bin").loc[[b[2] for b in bins]]
        x = np.arange(4) + (k - 1.5) * w; a.bar(x, s.sign_acc, w, color=c, label=models[m]); a.errorbar(x, s.sign_acc, yerr=[s.sign_acc - s.ci_lo, s.ci_hi - s.sign_acc], fmt="none", ecolor="k", elinewidth=0.5, capsize=1.5)
    a.axhline(0.5, color="k", lw=0.4); a.axhline(ref[cl]["zero"], color="#d62728", lw=0.6, ls="--", label=f"zero point, different RBPs ({ref[cl]['zero']:.2f})")
    a.axhline(ref[cl]["ceil"], color="#9467bd", lw=0.6, ls="--", label=f"same-RBP independent run ({ref[cl]['ceil']:.2f})")
    a.set_xticks(range(4)); a.set_xticklabels([f"{b[2]}\nn={int(S[(S.cell_line == cl) & (S.model == 'v3_A') & (S.repro_bin == b[2])].n_rbp.iloc[0])}" for b in bins]); a.set_ylim(0.4, 0.8)
    a.set_xlabel("RBP reproducibility bin"); a.set_ylabel("mean direction accuracy (95% CI)"); a.legend(frameon=False, fontsize=5.5, loc="upper left", ncol=2)
    a = ax[2, j]; n = N[N.cell_line == cl]; x = np.arange(len(n))
    a.bar(x - 0.3, n.v3_A, 0.2, color="#1f77b4", label="v3 full"); a.bar(x - 0.1, n.v3_none, 0.2, color="#2ca02c", label="event-only"); a.bar(x + 0.1, n.v3_shuf, 0.2, color="#7f7f7f", label="shuffled RBP"); a.bar(x + 0.3, n.agree_0_1, 0.2, color="#e7cfa0", label="own reproducibility")
    a.axhline(0.5, color="k", lw=0.4); a.set_xticks(x); a.set_xticklabels(n.target, rotation=60, ha="right", fontsize=6); a.set_ylim(0.1, 1.0); a.set_ylabel("direction accuracy"); a.legend(frameon=False, fontsize=5.5, ncol=2, loc="upper right")
    k = (n.group == "core spliceosome").sum(); a.axvline(k - 0.5, color="grey", lw=0.5, ls=":"); a.text(k / 2 - 0.5, 0.13, "core spliceosome", ha="center", fontsize=6); a.text(k + (len(n) - k) / 2 - 0.5, 0.13, "sequence-specific", ha="center", fontsize=6)
for i, lab in enumerate("ABC"): ax[i, 0].text(-0.14, 1.04, lab, transform=ax[i, 0].transAxes, fontsize=10, fontweight="bold")
plt.tight_layout(); plt.savefig(f"{M}/fig_stratified_v3.png"); plt.savefig(f"{M}/fig_stratified_v3.pdf")
pd.set_option("display.width", 250)
print(S[S.model.isin(["v3_A", "v3_none", "v3_shuf", "v2_A"])].pivot_table(index=["cell_line", "repro_bin"], columns="model", values="sign_acc").round(3).to_string())
print(); print(S[S.model == "v3_A"][["cell_line", "repro_bin", "n_rbp", "sign_acc", "ci_lo", "ci_hi", "diff_vs_shuffled", "diff_ci_lo", "diff_ci_hi"]].round(3).to_string(index=False))
print(); print(pd.DataFrame(ref).T.round(3).to_string()); print(); print(N.round(3).to_string(index=False))
print("files:", f"{M}/stratified_results.tsv", f"{M}/stratified_named_rbps.tsv", f"{M}/fig_stratified_v3.png", f"{M}/fig_stratified_v3.pdf")
