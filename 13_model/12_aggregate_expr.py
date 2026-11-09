"""Aggregate the gene-expression layer runs (13_model/runs_expr/<cell>_f<fold>_<tier><tag>/per_run.tsv) into one comparison:
tier A, tier B, gene-only control, shuffled-RBP control, and the training-run mean baseline, per RBP and per cell line, plus
paired bootstrap differences and stratification by the RBP's own expression reproducibility (expr_within_run_<cell>.tsv split-half r).
Outputs: 13_model/expr_vs_controls_by_rbp.tsv, expr_results_summary.tsv; stdout -> logs/expr_comparison.txt"""
import glob, numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; M = f"{PB}/13_model"; rng = np.random.default_rng(0); pd.set_option("display.width", 250)
V = {"A": ("_expr_s0", "A"), "B": ("_expr_s0", "B"), "none": ("_expr_ctrl_none_s0", "A"), "shuf": ("_expr_ctrl_shuf_s0", "A")}
def boot(x, n=2000):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 2: return np.nan, np.nan, np.nan
    b = rng.choice(x, (n, len(x))).mean(1); return x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5)
rows, summ = [], []
for cl in ("HepG2", "K562"):
    D = None
    for v, (tag, tier) in V.items():
        fs = sorted(glob.glob(f"{M}/runs_expr/{cl}_f[0-4]_{tier}{tag}/per_run.tsv"))
        if not fs: print(cl, v, "missing"); continue
        p = pd.concat([pd.read_csv(f, sep="\t") for f in fs]); by = p.groupby("target").agg(r_all=("r_all", "mean"), r_sig=("r_sig", "mean"), sign_sig=("sign_sig", "mean"), r_all_resid=("r_all_resid", "mean"), n_sig=("n_sig", "sum"),
                                                                                 r_all_meanbase=("r_all_meanbase", "mean"), r_sig_meanbase=("r_sig_meanbase", "mean"), sign_sig_meanbase=("sign_sig_meanbase", "mean"), mean_abs_pred=("mean_abs_pred", "mean"))
        by.columns = [f"{v}_{c}" if not c.endswith("meanbase") and c != "n_sig" else c for c in by.columns]; D = by if D is None else D.join(by[[c for c in by.columns if c.startswith(v + "_")]])
        print(cl, v, "folds", len(fs), "RBPs", len(by))
    w = pd.read_csv(f"{M}/expr_within_run_{cl}.tsv", sep="\t").groupby("target").r_split_half.mean(); D["r_split_half"] = w.reindex(D.index).values; D["cell_line"] = cl; D["target"] = D.index
    D["repro_bin"] = pd.cut(D.r_split_half, [-1, 0.2, 0.4, 0.6, 1.01], labels=["<0.2", "0.2-0.4", "0.4-0.6", ">0.6"]); rows.append(D)
    print(f"\n== {cl}: RBP-mean held-out metrics (expression layer) ==")
    cols = [c for c in D.columns if any(c.endswith(x) for x in ("_r_all", "_r_sig", "_sign_sig", "_r_all_resid"))] + ["r_all_meanbase", "r_sig_meanbase", "sign_sig_meanbase"]
    print(D[cols].mean().round(4).to_string())
    print(f"\n-- by RBP expression reproducibility (split-half r) --"); print(D.groupby("repro_bin", observed=True)[["A_r_all", "A_sign_sig", "none_r_all", "shuf_r_all", "r_all_meanbase", "A_r_all_resid"]].agg(["mean", "count"]).round(3).to_string())
    print("\n-- paired differences (RBP level, bootstrap 95% CI) --")
    for x, y in [("A_r_all", "shuf_r_all"), ("A_r_all", "none_r_all"), ("A_r_all", "r_all_meanbase"), ("A_r_all", "B_r_all"), ("A_sign_sig", "shuf_sign_sig"), ("A_sign_sig", "sign_sig_meanbase"), ("A_r_sig", "r_sig_meanbase")]:
        m, l, h = boot(D[x] - D[y]); print(f"{cl} {x} - {y}: {m:+.4f} [{l:+.4f},{h:+.4f}] n={int(np.isfinite(D[x] - D[y]).sum())}")
        summ.append(dict(cell_line=cl, a=x, b=y, diff=m, ci_lo=l, ci_hi=h))
    print("\n-- named RBPs --"); nm = ["SF3B1", "SF3A1", "U2AF2", "PRPF8", "PTBP1", "SRSF1", "HNRNPC", "UPF1", "XRN2", "EXOSC5", "LARP4", "PABPC4", "DDX6"]
    print(D[D.target.isin(nm)][["target", "r_split_half", "A_r_all", "B_r_all", "none_r_all", "shuf_r_all", "r_all_meanbase", "A_r_all_resid", "A_sign_sig", "shuf_sign_sig", "sign_sig_meanbase", "n_sig"]].round(3).to_string(index=False))
R = pd.concat(rows); R.to_csv(f"{M}/expr_vs_controls_by_rbp.tsv", sep="\t", index=False, float_format="%.4f"); pd.DataFrame(summ).to_csv(f"{M}/expr_results_summary.tsv", sep="\t", index=False, float_format="%.4f")
print("\nfiles:", f"{M}/expr_vs_controls_by_rbp.tsv", f"{M}/expr_results_summary.tsv")
