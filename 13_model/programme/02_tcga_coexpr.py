"""Observational prior: correlation between RBP expression and event PSI across TCGA primary tumours, within project, Fisher-z averaged over projects.
PSI: Snaptron junction counts already fetched for the dictionary events (14_inverse/tcga_psi, coverage inc+skip >= 10). Expression: UCSC Xena TOIL
log2(TPM+0.001). Two versions: raw, and adj = RBP expression residualised on the top 10 expression PCs of the project (5,000 most variable genes),
which absorbs purity / proliferation / subtype axes. No ENCODE data enters. Sign convention: stored value = r(PSI, expression); a knockdown is
expected to shift PSI by -r. Output: 06_features/coexpr/tcga_coexpr.npz {targets, events, raw, adj, n_projects}."""
import numpy as np, pandas as pd, glob, pyarrow.parquet as pq
PB = "/well/blagden/users/eic789/PerturbBridge"; O = f"{PB}/06_features/coexpr"
ids = pd.read_csv(f"{PB}/06_features/rbp_functional/rbp_id_map.tsv", sep="\t"); T = ids.target.values
sm = pd.read_csv(f"{PB}/14_inverse/tcga_samples.tsv", sep="\t", low_memory=False); sm = sm[sm.cgc_sample_sample_type.astype(str).str.startswith("Primary")].copy(); sm["xena"] = sm.tcga_barcode.str[:15]; sm = sm.drop_duplicates("xena")
print("reading expression", flush=True); X = pd.read_csv(f"{O}/tcga_RSEM_gene_tpm.gz", sep="\t", index_col=0).astype(np.float32); X.index = X.index.str.split(".").str[0]; X = X[~X.index.duplicated()]
pm = pd.read_csv(f"{O}/gencode.v23.annotation.gene.probemap", sep="\t"); pm["ensg"] = pm.id.str.split(".").str[0]; s2e = pm.drop_duplicates("gene").set_index("gene").ensg
gene = [e if isinstance(e, str) and e in X.index else s2e.get(t, s2e.get(s, None)) for t, e, s in zip(T, ids.ensg, ids.symbol)]; ok = np.array([g is not None and g in X.index for g in gene])
print("RBPs with expression", int(ok.sum()), "of", len(T), "| missing", list(T[~ok]), flush=True)
sm = sm[sm.xena.isin(X.columns)]; rail2col = dict(zip(sm.rail_id, range(len(sm)))); print("matched primary tumours", len(sm), flush=True)
parts = sorted(glob.glob(f"{PB}/14_inverse/tcga_psi/part_*.parquet")); ev = np.array(sorted(set().union(*[set(pq.read_table(p, columns=["event_id"]).column(0).unique().to_pylist()) for p in parts]))); epos = pd.Index(ev)
P = np.full((len(ev), len(sm)), np.nan, np.float32)
for p in parts:
    d = pq.read_table(p).to_pandas(); d = d[(d.inc + d.skip >= 10) & d.rail_id.isin(rail2col)]; P[epos.get_indexer(d.event_id), d.rail_id.map(rail2col).values] = d.psi.values; print(p.split("/")[-1], len(d), flush=True)
print("PSI matrix", P.shape, "measured %.3f" % np.isfinite(P).mean(), flush=True)
zs = {k: np.zeros((len(T), len(ev)), np.float64) for k in ("raw", "adj")}; ws = np.zeros(len(ev)); npj = np.zeros(len(ev), int)
for proj, g in sm.groupby("study"):
    if len(g) < 100: continue
    col = g.rail_id.map(rail2col).values; Xp = X[g.xena.values]; n = len(g)
    E = np.zeros((len(T), n), np.float32); E[ok] = Xp.loc[[x for x, o in zip(gene, ok) if o]].values
    hv = Xp[(Xp.mean(1) > 0)]; hv = hv.loc[hv.var(1).nlargest(5000).index].values; hv = hv - hv.mean(1, keepdims=True); Q = np.linalg.svd(hv, full_matrices=False)[2][:10].T
    def z(A): A = A - A.mean(1, keepdims=True); return A / (A.std(1, keepdims=True) + 1e-9)
    Ev = {"raw": z(E), "adj": z(z(E) - (z(E) @ Q) @ Q.T)}
    Pp = P[:, col]; M = np.isfinite(Pp).astype(np.float32); P0 = np.where(M > 0, Pp, 0).astype(np.float32); ne = M.sum(1); use = ne >= max(30, 0.3 * n)
    Sy = P0.sum(1); Syy = (P0 ** 2).sum(1); vy = Syy - Sy ** 2 / np.maximum(ne, 1); use &= vy > 1e-6
    for k, Ek in Ev.items():
        Sx = M @ Ek.T; Sxx = M @ (Ek ** 2).T; Sxy = P0 @ Ek.T; nn = np.maximum(ne, 1)[:, None]
        r = (Sxy - Sx * Sy[:, None] / nn) / np.sqrt(np.maximum(Sxx - Sx ** 2 / nn, 1e-12) * np.maximum(vy, 1e-12)[:, None]); r = np.clip(np.where(use[:, None], r, 0), -0.999, 0.999)
        zs[k] += (np.arctanh(r) * np.where(use, ne - 3, 0)[:, None]).T
    ws += np.where(use, ne - 3, 0); npj += use; print(proj, n, "events used", int(use.sum()), flush=True)
res = {k: np.tanh(v / np.maximum(ws, 1)[None]).astype(np.float32) for k, v in zs.items()}
for k in res: res[k][~ok] = np.nan
np.savez_compressed(f"{O}/tcga_coexpr.npz", targets=T, events=ev, n_projects=npj, has=ok, **res); print("events with >= 5 projects", int((npj >= 5).sum())); print("DONE")
