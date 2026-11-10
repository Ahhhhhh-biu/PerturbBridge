"""Per-programme paired bootstrap from the records of 03_programme_model.py: Pearson r of predicted vs measured score on each aligned programme,
model minus the mean of the permuted-descriptor controls, 2,000 resamples over held-out RBPs. A programme passes when r exceeds the 95th percentile
of the permuted controls and the bootstrap 95 % CI of the difference excludes 0. Usage: python 05_programme_report.py <tag>"""
import numpy as np, pandas as pd, sys
PB = "/well/blagden/users/eic789/PerturbBridge"; O = f"{PB}/13_model/programme"; tag = sys.argv[1] if len(sys.argv) > 1 else ""; R = pd.read_parquet(f"{O}/records{tag}.parquet"); K = sum(c.startswith("p") and c[1:].isdigit() for c in R.columns); rng = np.random.default_rng(1); rows = []
def rowr(A, B): A = A - A.mean(1, keepdims=True); B = B - B.mean(1, keepdims=True); return (A * B).sum(1) / np.sqrt((A ** 2).sum(1) * (B ** 2).sum(1) + 1e-30)
for (c, s, m), g in R[R.method != "-"].groupby(["cell", "set", "method"]):
    w = g[g.variant == "real"].set_index("target"); tg = w.index.values; pv = [x.set_index("target").reindex(tg) for v, x in g[g.variant != "real"].groupby("variant")]; idx = rng.integers(0, len(tg), (2000, len(tg)))
    for j in range(1, K + 1):
        t = w[f"t{j}"].values; r = np.corrcoef(w[f"p{j}"].values, t)[0, 1]; rp = [np.corrcoef(x[f"p{j}"].values, t)[0, 1] for x in pv]; d = rowr(w[f"p{j}"].values[idx], t[idx]) - np.mean([rowr(x[f"p{j}"].values[idx], t[idx]) for x in pv], 0)
        rows.append(dict(cell=c, set=s, method=m, programme=j, r=r, perm_p95=np.percentile(rp, 95), diff=d.mean(), lo=np.percentile(d, 2.5), hi=np.percentile(d, 97.5), passed=bool(r > np.percentile(rp, 95) and np.percentile(d, 2.5) > 0)))
T = pd.DataFrame(rows); T.to_csv(f"{O}/results_programme_bootstrap{tag}.tsv", sep="\t", index=False, float_format="%.3f")
n = T.groupby(["cell", "set", "method"]).apply(lambda g: pd.Series({"P1_passed": bool(g.passed[g.programme == 1].iloc[0]), "n_passed_P2_P10": int(g.passed[g.programme > 1].sum()), "which": ",".join(map(str, g.programme[g.passed & (g.programme > 1)]))})).reset_index()
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 200); print(n.to_string(index=False)); print("DONE report")
