"""Benchmark package summary. Compares, on identical held-out rows (clean-target rows of the RBP-held-out folds):
  v3 deep model (tier A), its shuffled-RBP and event-only controls, the ridge baseline (+ shuffled), and the v4 1D-CNN variant (+ shuffled).
Per-RBP metrics: direction accuracy on clean rows, Pearson r, Spearman rho, MAE on clean rows, MAE on null rows (calibration).
Paired RBP-level bootstrap differences (deep - ridge, deep - shuffled, CNN - deep). Difficulty analysis by RBP reproducibility
(within-run direction agreement) and by sequence distance to the training set (max BLASTP identity to a training-fold RBP).
Outputs: 13_model/benchmark/master_table.tsv, master_by_rbp.tsv, difficulty_table.tsv, fig_benchmark.png/.pdf; stdout -> logs/benchmark_metrics.txt
"""
import glob, os, numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

PB = "/well/blagden/users/eic789/PerturbBridge"
M = f"{PB}/13_model"
B = f"{M}/benchmark"
rng = np.random.default_rng(0)
pd.set_option("display.width", 260)
for f in ("Arial.ttf", "Arial-Bold.ttf"):
    font_manager.fontManager.addfont(f"/well/blagden/users/eic789/.fonts/{f}")
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6, "savefig.dpi": 300, "pdf.fonttype": 42})

MODELS = {
    "v3 deep": "A_v3_s0",
    "v3 shuffled RBP": "A_v3_ctrl_shuf_s0",
    "v3 event-only": "A_v3_ctrl_none_s0",
    "ridge": "A_ridge_s0",
    "ridge shuffled RBP": "A_ridge_ctrl_shuf_s0",
    "v4 CNN": "A_v4cnn_s0",
    "v4 CNN shuffled RBP": "A_v4cnn_ctrl_shuf_s0",
}
rep = pd.read_csv(f"{M}/ceiling_within_run_by_rbp.tsv", sep="\t")[["cell_line", "target", "agree_0_1"]]
near = pd.read_csv(f"{B}/nearest_training_rbp.tsv", sep="\t")[["target", "max_ident_train", "near_bin"]]


def per_rbp_metrics(cl, suffix):
    """Collect held-out predictions of one model over the five folds and compute per-RBP metrics."""
    z = np.load(f"{M}/tensors_{cl}.npz", allow_pickle=True)
    clean = np.load(f"{M}/clean_flags_{cl}.npz")["row_clean"]
    y = z["row_y"]
    run_of_row = z["row_r"]
    target_of_run = z["targets"][z["run_target"]]
    parts = []
    for fold in range(5):
        f = f"{M}/runs/{cl}_f{fold}_{suffix}/pred_heldout.npz"
        if not os.path.exists(f):
            return None
        p = np.load(f)
        parts.append(pd.DataFrame({"idx": p["idx"], "pred": p["pred"]}))
    d = pd.concat(parts, ignore_index=True)
    d["y"] = y[d.idx.values]
    d["clean"] = clean[d.idx.values]
    d["target"] = target_of_run[run_of_row[d.idx.values]]
    out = []
    for tg, g in d.groupby("target"):
        s = g[g.clean]
        n = g[~g.clean]
        if len(s) < 20:
            continue
        rho = stats.spearmanr(s.y, s.pred).correlation if s.pred.std() > 0 else np.nan
        r = np.corrcoef(s.y, s.pred)[0, 1] if s.pred.std() > 0 else np.nan
        out.append(dict(cell_line=cl, target=tg, n_clean=len(s), sign_acc=float((np.sign(s.y) == np.sign(s.pred)).mean()),
                        pearson=r, spearman=rho, mae_clean=float(np.abs(s.y - s.pred).mean()), mae_null=float(np.abs(n.y - n.pred).mean()),
                        mean_abs_pred_clean=float(np.abs(s.pred).mean())))
    return pd.DataFrame(out)


def boot(x, n=2000):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return np.nan, np.nan, np.nan
    b = rng.choice(x, (n, len(x))).mean(1)
    return x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5)


tables = {}
for name, suffix in MODELS.items():
    parts = [per_rbp_metrics(cl, suffix) for cl in ("HepG2", "K562")]
    if any(p is None for p in parts):
        print("missing runs for", name)
        continue
    t = pd.concat(parts, ignore_index=True)
    t["model"] = name
    tables[name] = t
R = pd.concat(tables.values(), ignore_index=True).merge(rep, on=["cell_line", "target"], how="left").merge(near, on="target", how="left")
R["repro_bin"] = pd.cut(R.agree_0_1, [0, 0.55, 0.6, 0.7, 1.01], labels=["<0.55", "0.55-0.6", "0.6-0.7", ">0.7"])
R.to_csv(f"{B}/master_by_rbp.tsv", sep="\t", index=False, float_format="%.4f")

# master table: RBP-mean metrics per model and cell line with bootstrap CI on direction accuracy
rows = []
for (model, cl), g in R.groupby(["model", "cell_line"]):
    m, lo, hi = boot(g.sign_acc)
    rows.append(dict(model=model, cell_line=cl, n_rbp=len(g), sign_acc=m, sign_ci_lo=lo, sign_ci_hi=hi, pearson=np.nanmean(g.pearson),
                     spearman=np.nanmean(g.spearman), mae_clean=g.mae_clean.mean(), mae_null=g.mae_null.mean(), mean_abs_pred_clean=g.mean_abs_pred_clean.mean()))
T = pd.DataFrame(rows)
T.to_csv(f"{B}/master_table.tsv", sep="\t", index=False, float_format="%.4f")
print("== master table (RBP mean on clean held-out rows) ==")
print(T.round(4).to_string(index=False))

# paired differences
print("\n== paired RBP-level differences (bootstrap 95% CI) ==")
W = R.pivot_table(index=["cell_line", "target"], columns="model", values="sign_acc")
Wr = R.pivot_table(index=["cell_line", "target"], columns="model", values="spearman")
pairs = [("v3 deep", "ridge"), ("v3 deep", "v3 shuffled RBP"), ("ridge", "ridge shuffled RBP"), ("v4 CNN", "v3 deep"), ("v4 CNN", "v4 CNN shuffled RBP")]
diff_rows = []
for a, b in pairs:
    if a not in W.columns or b not in W.columns:
        continue
    for cl in ("HepG2", "K562"):
        w = W.loc[cl]
        m, lo, hi = boot(w[a] - w[b])
        ms, los, his = boot(Wr.loc[cl][a] - Wr.loc[cl][b])
        diff_rows.append(dict(cell_line=cl, comparison=f"{a} - {b}", d_sign_acc=m, ci_lo=lo, ci_hi=hi, d_spearman=ms, sp_ci_lo=los, sp_ci_hi=his, n=int(np.isfinite(w[a] - w[b]).sum())))
D = pd.DataFrame(diff_rows)
D.to_csv(f"{B}/paired_differences.tsv", sep="\t", index=False, float_format="%.4f")
print(D.round(4).to_string(index=False))

# difficulty analysis
print("\n== difficulty: by RBP reproducibility and by sequence distance to the training set ==")
keep = [m for m in ("v3 deep", "ridge", "v4 CNN", "v3 event-only", "ridge shuffled RBP", "v3 shuffled RBP") if m in tables]
diff = []
for key in ("repro_bin", "near_bin"):
    g = R[R.model.isin(keep)].groupby(["cell_line", key, "model"], observed=True).agg(n=("target", "size"), sign_acc=("sign_acc", "mean"), spearman=("spearman", "mean"), mae_clean=("mae_clean", "mean")).reset_index()
    g["stratum"] = key
    diff.append(g)
    print(g.pivot_table(index=["cell_line", key], columns="model", values="sign_acc", observed=True).round(3).to_string())
    print(g.pivot_table(index=["cell_line", key], columns="model", values="n", observed=True).iloc[:, :1].to_string())
pd.concat(diff).to_csv(f"{B}/difficulty_table.tsv", sep="\t", index=False, float_format="%.4f")

# figure: A sign accuracy per model with CI; B difficulty by reproducibility; C difficulty by distance to training set; D MAE clean vs null (calibration)
fig, ax = plt.subplots(2, 2, figsize=(7.2, 6.2))
colors = {"v3 deep": "#1f77b4", "v3 shuffled RBP": "#7f7f7f", "v3 event-only": "#2ca02c", "ridge": "#ff7f0e", "ridge shuffled RBP": "#c7a26b", "v4 CNN": "#9467bd", "v4 CNN shuffled RBP": "#c5b0d5"}
a = ax[0, 0]
order = [m for m in MODELS if m in tables]
for j, cl in enumerate(("HepG2", "K562")):
    t = T[T.cell_line == cl].set_index("model").loc[order]
    x = np.arange(len(order)) + (j - 0.5) * 0.4
    a.bar(x, t.sign_acc, 0.38, color=[colors[m] for m in order], alpha=1 if j == 0 else 0.6, label=cl)
    a.errorbar(x, t.sign_acc, yerr=[t.sign_acc - t.sign_ci_lo, t.sign_ci_hi - t.sign_acc], fmt="none", ecolor="k", elinewidth=0.5, capsize=1.5)
a.axhline(0.5, color="k", lw=0.5)
a.set_xticks(range(len(order)))
a.set_xticklabels(order, rotation=35, ha="right", fontsize=6)
a.set_ylim(0.45, 0.6)
a.set_ylabel("direction accuracy, clean held-out rows (RBP mean, 95% CI)")
a.set_title("Same held-out RBPs (dark HepG2, light K562)", fontsize=8, loc="left")
for i, (key, lab) in enumerate((("repro_bin", "RBP reproducibility (within-run direction agreement)"), ("near_bin", "max BLASTP identity to a training RBP"))):
    a = ax[0, 1] if i == 0 else ax[1, 0]
    g = R[R.model.isin(keep)].groupby([key, "model"], observed=True).sign_acc.agg(["mean", "sem", "size"]).reset_index()
    cats = list(g[key].cat.categories) if hasattr(g[key], "cat") else sorted(g[key].unique())
    for k, m in enumerate(keep):
        s = g[g.model == m].set_index(key).reindex(cats)
        a.errorbar(np.arange(len(cats)) + (k - (len(keep) - 1) / 2) * 0.11, s["mean"], yerr=s["sem"], fmt="o", ms=3, color=colors[m], label=m, lw=0.8)
    a.axhline(0.5, color="k", lw=0.5)
    ns = g[g.model == keep[0]].set_index(key).reindex(cats)["size"]
    a.set_xticks(range(len(cats)))
    a.set_xticklabels([f"{c}\nn={int(n) if np.isfinite(n) else 0}" for c, n in zip(cats, ns)], fontsize=6)
    a.set_xlabel(lab)
    a.set_ylabel("direction accuracy (mean, s.e.m.)")
    a.legend(frameon=False, fontsize=6)
    a.set_title("Which RBPs are hard" if i == 0 else "Distance to the training set does not predict difficulty", fontsize=8, loc="left")
a = ax[1, 1]
for m in order:
    t = R[R.model == m]
    a.scatter(t.mae_null, t.mae_clean, s=5, lw=0, color=colors[m], label=m, alpha=0.7)
a.set_xlabel("MAE on null rows (should be small)")
a.set_ylabel("MAE on clean rows")
a.set_title("Calibration: error on null vs clean rows, per RBP", fontsize=8, loc="left")
a.legend(frameon=False, fontsize=5.5, ncol=2)
for i, (r_, c_) in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
    ax[r_, c_].text(-0.14, 1.06, "ABCD"[i], transform=ax[r_, c_].transAxes, fontsize=10, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{B}/fig_benchmark.png")
plt.savefig(f"{B}/fig_benchmark.pdf")
print("files:", f"{B}/master_table.tsv", f"{B}/master_by_rbp.tsv", f"{B}/paired_differences.tsv", f"{B}/difficulty_table.tsv", f"{B}/fig_benchmark.png")
