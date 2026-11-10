"""Programme-level forward model on held-out RBPs (design and decision rules: 13_model/programme/PLAN.md).
Per cell line and fold: programme basis = first K right singular vectors of the centred training matrix; target = scores of each RBP on the basis;
predictors = multi-output ridge and 5-nearest-neighbour mean from RBP descriptor blocks; controls = descriptors permuted across RBPs (NPERM times),
training-mean (zero) predictor, cross-cell-line oracle (the same RBP's measured scores in the other cell line).
Usage: python 03_programme_model.py [--k 10] [--nperm 20] [--tag ""] [--no_cx]
Outputs: 13_model/programme/results_summary<tag>.tsv, results_components<tag>.tsv, results_strata<tag>.tsv, per_rbp<tag>.tsv"""
import numpy as np, pandas as pd, os, argparse
from scipy.optimize import linear_sum_assignment
ap = argparse.ArgumentParser(); ap.add_argument("--k", type=int, default=10); ap.add_argument("--nperm", type=int, default=20); ap.add_argument("--tag", default=""); ap.add_argument("--nboot", type=int, default=2000); ap.add_argument("--no_cx", action="store_true", help="ignore the TCGA co-variation block"); a = ap.parse_args()
PB = "/well/blagden/users/eic789/PerturbBridge"; O = f"{PB}/13_model/programme"; K = a.k; rng = np.random.default_rng(20261110)
z = np.load(f"{O}/matrices.npz", allow_pickle=True); D = np.load(f"{PB}/06_features/rbp_functional/rbp_descriptors.npz", allow_pickle=True); ev = z["events"]
fold = pd.read_csv(f"{PB}/11_splits/rbp_folds.tsv", sep="\t").set_index("target").fold; dpos = {t: i for i, t in enumerate(D["targets"])}
rep = pd.read_csv(f"{PB}/14_inverse/dict_rbps.tsv", sep="\t").set_index(["cell_line", "target"]).repro_class
cxf = f"{PB}/06_features/coexpr/tcga_coexpr.npz"; CX = np.load(cxf, allow_pickle=True) if os.path.exists(cxf) and not a.no_cx else None
BLOCKS = ["seq", "ps", "net", "dep"] + (["cx"] if CX is not None else []); cum = {f"S{i}": BLOCKS[:i + 1] for i in range(len(BLOCKS))}
SETS = dict(cum); SETS.update({f"{b}_only": [b] for b in BLOCKS[1:]}); SETS.update({"seq_net": ["seq", "net"], "seq_net_dep": ["seq", "net", "dep"]}); LAMS = [0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300]
def P(*x): print(*x, flush=True)
def block(name, idx):
    X = D[name][idx].astype(np.float64); has = D["has_" + name][idx]; mu = X[has].mean(0); sd = X[has].std(0) + 1e-9; X = (X - mu) / sd; X[~has] = 0
    if not has.all(): X = np.hstack([X, (has - has.mean())[:, None] / (has.std() + 1e-9)])
    return X / np.sqrt(X.shape[1])
def ridge_fit(X, S, lam): mx = X.mean(0); ms = S.mean(0); Xc = X - mx; W = np.linalg.solve(Xc.T @ Xc + lam * np.eye(X.shape[1]), Xc.T @ (S - ms)); return lambda Z: (Z - mx) @ W + ms
def knn_pred(Xtr, Str, Xte, k=5):
    A = Xtr / (np.linalg.norm(Xtr, axis=1, keepdims=True) + 1e-12); B = Xte / (np.linalg.norm(Xte, axis=1, keepdims=True) + 1e-12); nn = np.argsort(-(B @ A.T), 1)[:, :k]; return Str[nn].mean(1)
def r2(S, Sh): return 1 - ((S - Sh) ** 2).sum() / (S ** 2).sum()
def rowr(A, B): A = A - A.mean(1, keepdims=True); B = B - B.mean(1, keepdims=True); return (A * B).sum(1) / np.sqrt((A ** 2).sum(1) * (B ** 2).sum(1) + 1e-30)
recs = []          # one record per (cell, set, method, variant, target)
for c, o in (("HepG2", "K562"), ("K562", "HepG2")):
    T = z[f"targets_{c}"]; Y = z[f"Y_{c}"].astype(np.float64); C = z[f"C_{c}"]; f = fold.reindex(T).values; di = np.array([dpos[t] for t in T]); n = len(T)
    To = z[f"targets_{o}"]; Yo = z[f"Y_{o}"].astype(np.float64); opos = {t: i for i, t in enumerate(To)}; fo = fold.reindex(To).values
    Vref = np.linalg.svd(np.where(np.isfinite(Y), Y, np.nanmean(Y, 0)) - np.nanmean(Y, 0), full_matrices=False)[2][:K]
    XB = {b: block(b, di) for b in BLOCKS if b != "cx"}
    if CX is not None:
        cpos = pd.Index(CX["events"]).get_indexer(ev); ctar = pd.Index(CX["targets"]).get_indexer(T); Cm = -CX["adj"][ctar][:, np.where(cpos >= 0, cpos, 0)].astype(np.float64); Cm[:, cpos < 0] = 0; hascx = np.isfinite(Cm).all(1); Cm[~np.isfinite(Cm)] = 0
    perms = [rng.permutation(n) for _ in range(a.nperm)]
    for k in range(5):
        te = f == k; tr = ~te; mu = np.nanmean(Y[tr], 0); Yc = np.where(np.isfinite(Y), Y, mu) - mu; V = np.linalg.svd(Yc[tr], full_matrices=False)[2][:K]; S = Yc @ V.T; sd = S[tr].std(0)
        A = np.abs(V @ Vref.T); ri, ci = linear_sum_assignment(-A); sg = np.sign((V * Vref[ci]).sum(1)); 
        def aligned(Sx): out = np.zeros_like(Sx); out[:, ci] = Sx[:, ri] * sg[ri]; return out
        if CX is not None:
            cf = (Cm - Cm[tr & hascx].mean(0)) @ V.T; cf = cf / (cf[tr & hascx].std(0) + 1e-9); cf[~hascx] = 0; XB["cx"] = np.hstack([cf, ((hascx - hascx.mean()) / (hascx.std() + 1e-9))[:, None]]) / np.sqrt(K + 1)
        # cross-cell-line oracle and reference reconstructions
        muo = np.nanmean(Yo[fo != k], 0); oi = np.array([opos.get(t, -1) for t in T]); So = np.full((n, K), np.nan); m = oi >= 0; So[m] = (np.where(np.isfinite(Yo[oi[m]]), Yo[oi[m]], muo) - muo) @ V.T
        mt = tr & m; slope = (So[mt] * S[mt]).sum(0) / (So[mt] ** 2).sum(0); preds = {("oracle_other_cell", "-", "real"): So * slope, ("zero", "-", "real"): np.zeros((n, K)), ("true_scores", "-", "real"): S}
        for sname, bl in SETS.items():
            X = np.hstack([XB[b] for b in bl])
            for var, Xv in [("real", X)] + [(f"perm{j}", X[p]) for j, p in enumerate(perms)]:
                Sz = S / sd; inner = [kk for kk in range(5) if kk != k]; score = []
                for lam in LAMS:
                    e = 0
                    for kk in inner: iv = f == kk; it = tr & ~iv; e += ((Sz[iv] - ridge_fit(Xv[it], Sz[it], lam)(Xv[iv])) ** 2).sum()
                    score.append(e)
                lam = LAMS[int(np.argmin(score))]; preds[(sname, "ridge", var)] = ridge_fit(Xv[tr], Sz[tr], lam)(Xv) * sd; preds[(sname, "knn", var)] = knn_pred(Xv[tr], S[tr], Xv)
        Yte = Y[te]; fin = np.isfinite(Yte); Cte = C[te] & fin; St = S[te]; Sal = aligned(St)
        for (sname, meth, var), Sh in preds.items():
            Sh = Sh[te]; rec = Sh @ V; ok = np.isfinite(Sh).all(1)
            er = np.array([np.corrcoef(rec[i, fin[i]], Yc[te][i, fin[i]])[0, 1] if ok[i] and rec[i].std() > 0 else np.nan for i in range(len(Sh))])
            hit = ((np.sign(mu + rec) == np.sign(Yte)) & Cte).sum(1); ncl = Cte.sum(1); Sha = aligned(np.where(np.isfinite(Sh), Sh, 0)); Sha[~ok] = np.nan
            for i, t in enumerate(T[te]):
                recs.append((c, sname, meth, var, t, k, ((St[i] - Sh[i]) ** 2).sum(), (St[i] ** 2).sum(), er[i], hit[i], ncl[i]) + tuple(Sha[i]) + tuple(Sal[i]))
        P(c, "fold", k, "done")
cols = ["cell", "set", "method", "variant", "target", "fold", "sse", "sst", "event_r", "hit", "n_clean"] + [f"p{j + 1}" for j in range(K)] + [f"t{j + 1}" for j in range(K)]
R = pd.DataFrame(recs, columns=cols); R["repro"] = [rep.get((c, t), "na") for c, t in zip(R.cell, R.target)]
R[R.variant == "real"].to_csv(f"{O}/per_rbp{a.tag}.tsv", sep="\t", index=False); R.to_parquet(f"{O}/records{a.tag}.parquet")
def metrics(g):
    g = g.dropna(subset=["p1"]); d = {"n": len(g), "r2_sub": 1 - g.sse.sum() / g.sst.sum(), "event_r": g.event_r.mean(), "dir_acc": (g.hit / g.n_clean)[g.n_clean >= 20].mean()}
    for j in range(K): d[f"r{j + 1}"] = np.corrcoef(g[f"p{j + 1}"], g[f"t{j + 1}"])[0, 1] if g[f"p{j + 1}"].std() > 0 else np.nan
    keep_ = g.t1.abs() <= g.t1.abs().quantile(0.95); d["rho1"] = g.p1.corr(g.t1, method="spearman") if g.p1.std() > 0 else np.nan; d["r1_trim"] = g.p1[keep_].corr(g.t1[keep_]) if g.p1.std() > 0 else np.nan
    return pd.Series(d)
M = R.groupby(["cell", "set", "method", "variant"]).apply(metrics).reset_index(); real = M[M.variant == "real"].drop(columns="variant"); perm = M[M.variant != "real"]
pm = perm.drop(columns="variant").groupby(["cell", "set", "method"]).agg(["mean", lambda x: np.nanpercentile(x, 95)]); pm.columns = [f"{a_}_perm_{'mean' if b == 'mean' else 'p95'}" for a_, b in pm.columns]
summ = real.merge(pm.reset_index(), on=["cell", "set", "method"], how="left")
# paired bootstrap over RBPs: model minus mean of permuted controls, for r1 and r2_sub; and consecutive cumulative sets
def boot(c, s, meth):
    g = R[(R.cell == c) & (R.set == s) & (R.method == meth)]; w = g[g.variant == "real"].set_index("target"); tg = w.index.values; pv = [g[g.variant == f"perm{j}"].set_index("target").reindex(tg) for j in range(a.nperm)]
    idx = rng.integers(0, len(tg), (a.nboot, len(tg))); t1 = w.t1.values[idx]
    def br(p): return rowr(p[idx], t1)
    def b2(x): return 1 - x.sse.values[idx].sum(1) / x.sst.values[idx].sum(1)
    d1 = br(w.p1.values) - np.mean([br(x.p1.values) for x in pv], 0); d2 = b2(w) - np.mean([b2(x) for x in pv], 0)
    return dict(r1_diff=d1.mean(), r1_lo=np.percentile(d1, 2.5), r1_hi=np.percentile(d1, 97.5), r2_diff=d2.mean(), r2_lo=np.percentile(d2, 2.5), r2_hi=np.percentile(d2, 97.5))
B = pd.DataFrame([dict(cell=c, set=s, method=m_, **boot(c, s, m_)) for c in ("HepG2", "K562") for s in SETS for m_ in ("ridge", "knn")]); summ = summ.merge(B, on=["cell", "set", "method"], how="left")
keep = ["cell", "set", "method", "n", "r1", "r1_perm_mean", "r1_perm_p95", "r1_diff", "r1_lo", "r1_hi", "rho1", "rho1_perm_p95", "r1_trim", "r1_trim_perm_p95", "r2_sub", "r2_sub_perm_mean", "r2_sub_perm_p95", "r2_diff", "r2_lo", "r2_hi", "event_r", "event_r_perm_mean", "dir_acc", "dir_acc_perm_mean"]
summ[keep].to_csv(f"{O}/results_summary{a.tag}.tsv", sep="\t", index=False, float_format="%.4f")
comp = summ[["cell", "set", "method"] + [x for j in range(K) for x in (f"r{j + 1}", f"r{j + 1}_perm_p95")]]; comp.to_csv(f"{O}/results_components{a.tag}.tsv", sep="\t", index=False, float_format="%.3f")
st = R[R.variant == "real"].groupby(["cell", "set", "method", "repro"]).apply(metrics).reset_index()[["cell", "set", "method", "repro", "n", "r1", "r2_sub", "event_r", "dir_acc"]]; st.to_csv(f"{O}/results_strata{a.tag}.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 200)
P("\n== summary (held-out RBPs pooled over folds)"); P(summ[keep].round(3).to_string(index=False))
P("\n== per-programme r (real | 95th percentile of permuted)"); P(comp.round(2).to_string(index=False))
P("\n== cumulative sets, paired bootstrap of subspace R2 (ridge)")
for c in ("HepG2", "K562"):
    names = list(cum)
    for s0, s1 in zip(names[:-1], names[1:]):
        w0 = R[(R.cell == c) & (R.set == s0) & (R.method == "ridge") & (R.variant == "real")].set_index("target"); w1 = R[(R.cell == c) & (R.set == s1) & (R.method == "ridge") & (R.variant == "real")].set_index("target").reindex(w0.index)
        idx = rng.integers(0, len(w0), (a.nboot, len(w0))); d = (1 - w1.sse.values[idx].sum(1) / w1.sst.values[idx].sum(1)) - (1 - w0.sse.values[idx].sum(1) / w0.sst.values[idx].sum(1))
        P(f"  {c} {s1} - {s0}: {d.mean():+.4f} [{np.percentile(d, 2.5):+.4f}, {np.percentile(d, 97.5):+.4f}]")
P("\n== strata (ridge, real)"); P(st[st.method.isin(["ridge", "-"])].round(3).to_string(index=False)); P("DONE")
