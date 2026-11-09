"""Baselines for the forward task on RBP-held-out folds (design 4.2), one cell line per call.
Per held-out run the prediction is a dPSI vector over the cell line's training events:
  zero        : 0
  event_mean  : mean dPSI of the event over training runs
  esm2_nn     : mean dPSI of the 3 training RBPs with highest ESM2 cosine to the held-out RBP
  genept_nn   : same with GenePT ada cosine
  eclip_nn    : (tier-A runs only) mean dPSI of the 3 training RBPs (with eCLIP in this cell line) whose eCLIP event profile
                (max signal over up/seg/dn) correlates best with the held-out RBP's profile
  ridge       : per-event ridge of dPSI on [ESM2 PCA64, GenePT PCA64, assay]; lambda chosen on an inner fold
Metrics per run: pearson over all covered events (r_all), pearson over the run's significant events (r_sig, FDR<0.05 & |dPSI|>0.1),
sign accuracy on significant events (sign_sig), MSE; sign accuracy by |dPSI| bin. Aggregation: mean per RBP, then mean over RBPs
with a bootstrap CI over RBPs; reported by tier.
Usage: python 02_baselines.py <cell_line> [dpsi|dpsi_cc]  -> 12_baselines/baselines_<cell_line>_<label>_per_run.tsv, _summary.tsv
(dpsi_cc = control-centred dPSI from 03_control_centering.py)"""
import sys, os, pickle, numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; U = f"{PB}/10_event_universe/encode_full"; O = f"{PB}/12_baselines"; cl = sys.argv[1]; LAB = sys.argv[2] if len(sys.argv) > 2 else "dpsi"; rng = np.random.default_rng(0)
tiers = pd.read_csv(f"{PB}/11_splits/run_tiers.tsv", sep="\t").set_index("run")
t = pd.read_parquet(f"{U}/train_table_{cl}.parquet", columns=["event_id", "run", "target", "dpsi", "fdr", LAB] if LAB != "dpsi" else ["event_id", "run", "target", "dpsi", "fdr"])
M = t.pivot(index="event_id", columns="run", values=LAB)
if os.environ.get("PB_TEST"): M = M.iloc[:3000]; O = f"{PB}/envs/tmp"
SIG = ((t.fdr < 0.05) & (t.dpsi.abs() > 0.1)).groupby([t.event_id, t.run]).first().unstack().reindex(index=M.index, columns=M.columns).fillna(False).astype(bool)
Y = M.values; runs = M.columns.to_list(); ev = M.index.to_list(); info = tiers.loc[runs]; tgt = info.target.values; fold = info.fold.values
print(cl, "events", len(ev), "runs", len(runs), file=sys.stderr, flush=True)
# RBP features
ix = pd.read_csv(f"{PB}/06_features/rbp_annotation/rbp_esm2_650M_index.tsv", sep="\t"); E = np.load(f"{PB}/06_features/rbp_annotation/rbp_esm2_650M_mean.npy"); esm = dict(zip(ix.target, E))
gp = pickle.load(open(f"{PB}/06_features/rbp_annotation/genept/GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle", "rb"))
targets = sorted(set(tgt)); miss_gp = [x for x in targets if x not in gp]; print("targets without GenePT:", miss_gp, file=sys.stderr)
def unit(X): return X / np.linalg.norm(X, axis=1, keepdims=True)
Fe = unit(np.array([esm[x] for x in targets])); Fg = unit(np.array([np.asarray(gp[x]) if x in gp else np.zeros(1536) for x in targets])); tpos = {x: i for i, x in enumerate(targets)}
def pca(X, k):
    Xc = X - X.mean(0); u, s, vt = np.linalg.svd(Xc, full_matrices=False); Z = Xc @ vt[:k].T; return Z / Z.std(0)
X = np.hstack([pca(Fe, 64), pca(Fg, 64)]); X = X[[tpos[x] for x in tgt]]; X = np.hstack([X, (info.assay.str.startswith("CRISPR").values[:, None]).astype(float)])
# eCLIP profiles (max signal over regions) for RBPs with eCLIP in this cell line
ec = pd.read_parquet(f"{U}/eclip_event_features.parquet"); ec = ec[(ec.cell_line == cl) & ec.event_id.isin(ev)]
ec["s"] = ec[["sig_up", "sig_seg", "sig_dn"]].max(1); P = ec.pivot_table(index="event_id", columns="rbp", values="s", aggfunc="max").reindex(ev).fillna(0)
Pz = (P - P.mean()) / (P.std() + 1e-9); eclip_rbps = set(P.columns)
def pearson_cols(A, B):  # column-wise pearson between matrices with matching rows, NaN-aware for A
    out = np.full(A.shape[1], np.nan)
    for j in range(A.shape[1]):
        m = np.isfinite(A[:, j]) & np.isfinite(B[:, j]); out[j] = np.corrcoef(A[m, j], B[m, j])[0, 1] if m.sum() > 2 else np.nan
    return out
def nn_pred(sim, train_cols, k=3):  # sim: similarity to each training target (dict target -> sim)
    top = sorted(sim, key=sim.get, reverse=True)[:k]; cols = [j for j in train_cols if tgt[j] in top]; return np.nanmean(Y[:, cols], axis=1)
rows = []; bins = [(0.05, 0.1), (0.1, 0.2), (0.2, 1.01)]
for k in range(5):
    te = np.where(fold == k)[0]; tr = np.where(fold != k)[0]; inner = tr[fold[tr] == (k + 1) % 5]; trin = tr[fold[tr] != (k + 1) % 5]
    mean_tr = np.nanmean(Y[:, tr], axis=1); tr_targets = sorted(set(tgt[tr]))
    # ridge: fit on imputed, event-centred Y; choose lambda on inner fold
    def ridge_fit(cols, lam):
        Yc = Y[:, cols].copy(); mu = np.nanmean(Yc, axis=1, keepdims=True); Yc = np.where(np.isfinite(Yc), Yc, mu) - mu
        Xc = X[cols]; xm = Xc.mean(0); Xc = Xc - xm; W = np.linalg.solve(Xc.T @ Xc + lam * np.eye(Xc.shape[1]), Xc.T @ Yc.T); return W, xm, mu[:, 0]
    best = None
    for lam in (10.0, 100.0, 1000.0, 10000.0):
        W, xm, mu = ridge_fit(trin, lam); pred = ((X[inner] - xm) @ W).T + mu[:, None]; sc = np.nanmean(pearson_cols(Y[:, inner], pred))
        if best is None or sc > best[0]: best = (sc, lam)
    W, xm, mu = ridge_fit(tr, best[1]); print(f"fold {k}: ridge lambda {best[1]} (inner r {best[0]:.3f})", file=sys.stderr, flush=True)
    for j in te:
        y = Y[:, j]; obs = np.isfinite(y); sig = SIG.values[:, j] & obs; tg = tgt[j]; tier = info.tier.values[j]
        preds = {"zero": np.zeros(len(ev)), "event_mean": mean_tr,
                 "esm2_nn": nn_pred({x: float(Fe[tpos[tg]] @ Fe[tpos[x]]) for x in tr_targets}, tr),
                 "genept_nn": nn_pred({x: float(Fg[tpos[tg]] @ Fg[tpos[x]]) for x in tr_targets}, tr),
                 "ridge": (X[j] - xm) @ W + mu}
        if tier == "A" and tg in eclip_rbps:
            cand = [x for x in tr_targets if x in eclip_rbps]
            if cand: preds["eclip_nn"] = nn_pred({x: float(np.corrcoef(Pz[tg], Pz[x])[0, 1]) for x in cand}, tr)
        for name, p in preds.items():
            p = np.where(np.isfinite(p), p, 0.0)  # NN predictions are NaN where no neighbour run covers the event: predict no change there
            d = dict(run=runs[j], target=tg, tier=tier, fold=k, method=name, n_obs=int(obs.sum()), n_sig=int(sig.sum()))
            d["r_all"] = np.corrcoef(y[obs], p[obs])[0, 1] if p[obs].std() > 0 else np.nan
            d["r_sig"] = np.corrcoef(y[sig], p[sig])[0, 1] if sig.sum() > 2 and p[sig].std() > 0 else np.nan
            d["sign_sig"] = (np.sign(y[sig]) == np.sign(p[sig])).mean() if sig.sum() else np.nan
            d["mse_all"] = float(np.mean((y[obs] - p[obs]) ** 2)); d["mse_sig"] = float(np.mean((y[sig] - p[sig]) ** 2)) if sig.sum() else np.nan
            for lo, hi in bins:
                b = obs & (np.abs(y) >= lo) & (np.abs(y) < hi); d[f"sign_{lo}_{hi}"] = (np.sign(y[b]) == np.sign(p[b])).mean() if b.sum() else np.nan
            rows.append(d)
R = pd.DataFrame(rows); R.to_csv(f"{O}/baselines_{cl}_{LAB}_per_run.tsv", sep="\t", index=False)
met = [c for c in R.columns if c.startswith(("r_", "sign_", "mse_"))]
byrbp = R.groupby(["method", "tier", "target"])[met].mean().reset_index()
def boot(g):
    out = {}
    for m in met:
        v = g[m].dropna().values
        if len(v) == 0: out[m] = "NA"; continue
        bs = [rng.choice(v, len(v)).mean() for _ in range(1000)]; out[m] = f"{v.mean():.3f} [{np.percentile(bs, 2.5):.3f}, {np.percentile(bs, 97.5):.3f}]"
    out["n_rbp"] = len(g); return pd.Series(out)
S = pd.concat([byrbp.groupby(["method", "tier"]).apply(boot), byrbp.assign(tier="all").groupby(["method", "tier"]).apply(boot)]).reset_index()
S.to_csv(f"{O}/baselines_{cl}_{LAB}_summary.tsv", sep="\t", index=False); print(S.to_string(index=False))
