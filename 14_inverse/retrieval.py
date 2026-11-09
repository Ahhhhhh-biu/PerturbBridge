"""Shared inverse-mapping helpers: dictionary loading, event-key construction from rMATS tables, per-sample PSI parsing, scoring and nulls."""
import numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; O = f"{PB}/14_inverse"; U = f"{PB}/10_event_universe/encode_full"
COORD = {"SE": ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
         "A3SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
         "A5SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
         "MXE": ["1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
         "RI": ["riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"]}
def load_dictionary(cells=("HepG2", "K562"), min_runs=1):
    """events x RBPs matrix of mean control-centred dPSI; RBPs present in both cell lines are averaged over cell lines."""
    mats = []
    for cl in cells:
        g = pd.read_parquet(f"{O}/dict_splicing_rbp_{cl}.parquet"); mats.append(g.pivot(index="event_id", columns="target", values="dpsi_cc"))
    if len(mats) == 1: return mats[0]
    idx = mats[0].index.union(mats[1].index); cols = mats[0].columns.union(mats[1].columns)
    A = mats[0].reindex(index=idx, columns=cols); B = mats[1].reindex(index=idx, columns=cols)
    return pd.DataFrame(np.nanmean(np.stack([A.values, B.values]), 0), index=idx, columns=cols)
def rbp_info():
    d = pd.read_csv(f"{O}/dict_rbps.tsv", sep="\t"); best = d.sort_values("agree_0_1", ascending=False).drop_duplicates("target").set_index("target")
    return best[["agree_0_1", "repro_class", "n_clean_events"]]
def rmats_psi_table(d, ev_key_to_id=None):
    """Read SE/A3SS/A5SS/MXE/RI JC tables of one rMATS run; return long table event_id x sample with psi and coverage (IJC+SJC)
    using rMATS's own per-sample IncLevel; samples are 'g1_<i>' / 'g2_<i>' in b1/b2 order."""
    if ev_key_to_id is None:
        ev = pd.read_parquet(f"{U}/events.parquet", columns=["key", "event_id"]); ev_key_to_id = dict(zip(ev.key, ev.event_id))
    out = []
    for t, cols in COORD.items():
        x = pd.read_csv(f"{d}/{t}.MATS.JC.txt", sep="\t", dtype=str)
        key = t + "|" + x["chr"] + "|" + x["strand"] + "|" + x[cols].astype(str).agg("|".join, axis=1); x["event_id"] = key.map(ev_key_to_id); x = x[x.event_id.notna()]
        for grp, inc, sk, lv in (("g1", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1"), ("g2", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")):
            I = x[inc].str.split(",", expand=True).astype(float); S = x[sk].str.split(",", expand=True).astype(float); L = x[lv].str.split(",", expand=True).replace("NA", np.nan).astype(float)
            cov = I + S; n = I.shape[1]
            for j in range(n):
                out.append(pd.DataFrame({"event_id": x.event_id.values, "sample": f"{grp}_{j}", "psi": L[j].values, "cov": cov[j].values}))
    return pd.concat(out, ignore_index=True)
def score_against(q, P, min_shared=50):
    """Pearson r of query vector q (pd.Series event->delta) with each column of P over pairwise-complete events. Returns (scores, n_shared)."""
    M = P.loc[P.index.intersection(q.index)]; qv = q.loc[M.index].values.astype(float); X = M.values
    ok = np.isfinite(X) & np.isfinite(qv)[:, None]; n = ok.sum(0); qx = np.where(ok, qv[:, None], 0.0); Xx = np.where(ok, X, 0.0)
    mq = qx.sum(0) / np.maximum(n, 1); mx = Xx.sum(0) / np.maximum(n, 1)
    cov = ((qx - mq) * (Xx - mx) * ok).sum(0); vq = (((qx - mq) ** 2) * ok).sum(0); vx = (((Xx - mx) ** 2) * ok).sum(0)
    r = cov / np.sqrt(np.maximum(vq * vx, 1e-12)); r[n < min_shared] = np.nan
    return pd.Series(r, index=M.columns), pd.Series(n, index=M.columns)
def ols_group_effect(Y, X):
    """Y: samples x events (NaN allowed), X: samples x p design with the group indicator in column 1 (column 0 = intercept).
    Returns per-event coefficient, t statistic and two-sided p of the group column (events with < p+3 complete samples -> NaN).
    Vectorised over events: per-event masked normal equations (XtX_e = X' diag(m_e) X) solved in a batch; identical to a per-event OLS."""
    from scipy import stats
    n_s, n_e = Y.shape; p = X.shape[1]; M = np.isfinite(Y); Yz = np.where(M, Y, 0.0); Mf = M.astype(float); nobs = Mf.sum(0)
    XtX = np.einsum("si,sj,se->eij", X, X, Mf); Xty = np.einsum("si,se->ei", X, Yz)
    ok = nobs >= p + 3; beta = np.full(n_e, np.nan); tval = np.full(n_e, np.nan); pval = np.full(n_e, np.nan)
    if ok.any():
        A = XtX[ok]; A = A + np.eye(p)[None] * 1e-10
        try: inv = np.linalg.inv(A)
        except np.linalg.LinAlgError: inv = np.linalg.pinv(A)
        B = np.einsum("eij,ej->ei", inv, Xty[ok])                                  # e x p coefficients
        res = (Yz[:, ok] - X @ B.T) * Mf[:, ok]; dof = nobs[ok] - p; s2 = (res ** 2).sum(0) / dof
        se = np.sqrt(s2 * inv[:, 1, 1]); good = s2 > 1e-20   # zero residual variance (all observed PSI identical): no test, t and p -> NaN
        beta[ok] = B[:, 1]; tt = np.where(good, B[:, 1] / np.where(good, se, 1.0), np.nan); tval[ok] = tt; pval[ok] = np.where(good, 2 * stats.t.sf(np.abs(tt), dof), np.nan)
    return beta, tval, pval
def ols_group_effect_slow(Y, X):
    """Reference per-event implementation (kept for checks)."""

    from scipy import stats
    n_s, n_e = Y.shape; p = X.shape[1]; beta = np.full(n_e, np.nan); tval = np.full(n_e, np.nan); pval = np.full(n_e, np.nan)
    ok = np.isfinite(Y); full = ok.all(0)
    # fast path for events with all samples
    if full.any():
        XtX_inv = np.linalg.pinv(X.T @ X); B = XtX_inv @ X.T @ Y[:, full]; res = Y[:, full] - X @ B; dof = n_s - p
        s2 = (res ** 2).sum(0) / dof; se = np.sqrt(s2 * XtX_inv[1, 1]); beta[full] = B[1]; tval[full] = B[1] / se; pval[full] = 2 * stats.t.sf(np.abs(tval[full]), dof)
    for j in np.where(~full)[0]:
        m = ok[:, j]
        if m.sum() < p + 3: continue
        Xi = X[m]; yi = Y[m, j]; XtX_inv = np.linalg.pinv(Xi.T @ Xi); b = XtX_inv @ Xi.T @ yi; res = yi - Xi @ b; dof = m.sum() - p
        s2 = (res ** 2).sum() / dof; se = np.sqrt(s2 * XtX_inv[1, 1]); beta[j] = b[1]; tval[j] = b[1] / se if se > 0 else np.nan; pval[j] = 2 * stats.t.sf(abs(tval[j]), dof) if se > 0 else np.nan
    return beta, tval, pval
def rank_report(scores, info, label, top=15):
    s = scores.dropna().sort_values(ascending=False); d = pd.DataFrame({"rbp": s.index, "score": s.values, "rank": np.arange(1, len(s) + 1)}).set_index("rbp").join(info)
    d["query"] = label; return d.reset_index()
def welch_nan(A, B):
    """Column-wise Welch t-test between rows of A and rows of B (NaN allowed). Returns (t, p, mean difference A - B).
    Same result as scipy.stats.ttest_ind(..., equal_var=False, nan_policy='omit') but vectorised over columns with masked sums."""
    from scipy import stats
    def msd(M):
        ok = np.isfinite(M); n = ok.sum(0).astype(float); Z = np.where(ok, M, 0.0); m = Z.sum(0) / np.maximum(n, 1)
        v = ((np.where(ok, M - m, 0.0)) ** 2).sum(0) / np.maximum(n - 1, 1); return n, m, v
    na, ma, va = msd(A); nb, mb, vb = msd(B); sa = va / np.maximum(na, 1); sb = vb / np.maximum(nb, 1); se = np.sqrt(sa + sb)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (ma - mb) / se; df = (sa + sb) ** 2 / (sa ** 2 / np.maximum(na - 1, 1) + sb ** 2 / np.maximum(nb - 1, 1))
    bad = (na < 2) | (nb < 2) | ~np.isfinite(t); t = np.where(bad, np.nan, t); p = np.where(bad, np.nan, 2 * stats.t.sf(np.abs(t), df))
    return t, p, ma - mb
def nan_totals(Y):
    """Per-column (count, sum, sum of squares) over finite entries of Y; feeds welch_split for repeated label permutations."""
    ok = np.isfinite(Y); Z = np.where(ok, Y, 0.0); return ok.sum(0).astype(float), Z.sum(0), (Z ** 2).sum(0)
def welch_split(Y, mask, tot):
    """Welch t-test of rows Y[mask] against rows Y[~mask] using precomputed totals (nan_totals(Y)); only the case rows are touched,
    so a permutation costs O(n_case x events). Returns (t, p, mean difference case - control). Agrees with welch_nan."""
    from scipy import stats
    A = Y[mask]; ok = np.isfinite(A); Z = np.where(ok, A, 0.0); na = ok.sum(0).astype(float); s1a = Z.sum(0); s2a = (Z ** 2).sum(0)
    n, s1, s2 = tot; nb = n - na; s1b = s1 - s1a; s2b = s2 - s2a
    with np.errstate(divide="ignore", invalid="ignore"):
        ma = s1a / na; mb = s1b / nb; va = np.maximum(s2a - na * ma ** 2, 0) / (na - 1); vb = np.maximum(s2b - nb * mb ** 2, 0) / (nb - 1)
        sa = va / na; sb = vb / nb; t = (ma - mb) / np.sqrt(sa + sb); df = (sa + sb) ** 2 / (sa ** 2 / (na - 1) + sb ** 2 / (nb - 1))
    bad = (na < 2) | (nb < 2) | ~np.isfinite(t); t = np.where(bad, np.nan, t); p = np.where(bad, np.nan, 2 * stats.t.sf(np.abs(t), df))
    return t, p, ma - mb
