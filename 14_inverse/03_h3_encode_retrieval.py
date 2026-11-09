"""Block 2a: H3 qualification inside ENCODE. Can a knockdown profile retrieve its own RBP from the dictionary?
Query = one run's clean events (replicate-consistent, control-independent) with their control-centred dPSI. Candidate score = Pearson r
between the query values and the candidate RBP's mean dpsi_cc over the shared events (>= MINSHARED events). Rank of the true RBP
among all candidate RBPs of the dictionary cell line (1 = best). Two settings:
  same   : same cell line, leave-one-run-out (true RBP keeps its other runs) -> only RBPs with >= 2 runs in that cell line
  cross  : query run from one cell line, dictionary = the other cell line (all runs) -> every RBP present in both cell lines
Null: ranks are uniform on 1..N under no information (expected top-5 rate 5/N); also an empirical null from 200 label permutations
of the dictionary columns per query. Family-level hit: true RBP's paralog group (11_splits/paralog_groups.tsv) in top-5.
Outputs: 14_inverse/h3_encode_retrieval.tsv (per query run), h3_encode_summary.tsv; stdout -> logs/h3_encode.txt"""
import numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; O = f"{PB}/14_inverse"; MINSHARED = 50; rng = np.random.default_rng(0); pd.set_option("display.width", 250)
rb = pd.read_csv(f"{O}/dict_rbps.tsv", sep="\t"); par = pd.read_csv(f"{PB}/11_splits/paralog_groups.tsv", sep="\t")
fam = {m: g for g, ms in zip(par.paralog_group, par.members) for m in str(ms).split(",")}   # member -> paralog group
runs, prof = {}, {}
for cl in ("HepG2", "K562"):
    t = pd.read_parquet(f"{O}/dict_splicing_runs_{cl}.parquet"); runs[cl] = t
    g = pd.read_parquet(f"{O}/dict_splicing_rbp_{cl}.parquet"); prof[cl] = g.pivot(index="event_id", columns="target", values="dpsi_cc")   # events x RBPs (NaN where not covered)
    print(cl, "profile matrix", prof[cl].shape)
def corr_cols(q, M):
    """Pearson r between query vector q (events,) and every column of M (events x RBPs) over pairwise-complete rows."""
    X = M.values; ok = np.isfinite(X) & np.isfinite(q)[:, None]; n = ok.sum(0)
    qx = np.where(ok, q[:, None], 0.0); Xx = np.where(ok, X, 0.0)
    mq = qx.sum(0) / np.maximum(n, 1); mx = Xx.sum(0) / np.maximum(n, 1)
    cov = ((qx - mq) * (Xx - mx) * ok).sum(0); vq = (((qx - mq) ** 2) * ok).sum(0); vx = (((Xx - mx) ** 2) * ok).sum(0)
    r = cov / np.sqrt(np.maximum(vq * vx, 1e-12)); r[n < MINSHARED] = np.nan; return pd.Series(r, index=M.columns), pd.Series(n, index=M.columns)
out = []
def evaluate(setting, qcl, dcl):
    t = runs[qcl]; P = prof[dcl]; info = rb[rb.cell_line == qcl].set_index("target")
    rt = t[["run", "target", "assay"]].drop_duplicates()
    if setting == "same": multi = rt.groupby("target").run.nunique(); rt = rt[rt.target.isin(multi[multi >= 2].index)]
    else: rt = rt[rt.target.isin(P.columns)]
    tq = t[t.clean]
    for r in rt.itertuples():
        q = tq[tq.run == r.run].set_index("event_id").dpsi_cc
        if len(q) < MINSHARED: continue
        M = P.loc[P.index.intersection(q.index)]; qv = q.loc[M.index].values.astype(float)
        if setting == "same":   # rebuild the true RBP's column without the query run
            oth = t[(t.target == r.target) & (t.run != r.run) & t.event_id.isin(M.index)].groupby("event_id").dpsi_cc.mean(); M = M.copy(); M[r.target] = np.nan; M.loc[oth.index, r.target] = oth.values
        s, n = corr_cols(qv, M); s = s.dropna()
        if r.target not in s.index: continue
        rank = int((s > s[r.target]).sum() + 1); N = len(s); top = s.sort_values(ascending=False)
        famhit = any(fam.get(x, x) == fam.get(r.target, r.target) for x in top.index[:5])
        # empirical null: permute column labels (equivalently the true RBP gets a random other RBP's score)
        null_rank = (rng.random(200) * N).astype(int) + 1
        out.append(dict(setting=setting, query_cell=qcl, dict_cell=dcl, run=r.run, target=r.target, assay=r.assay, n_query_events=len(q), n_candidates=N, rank=rank, rank_frac=rank / N,
                        score_true=s[r.target], score_top1=top.iloc[0], top1=top.index[0], top5=";".join(top.index[:5]), family_top5=famhit,
                        agree_0_1=info.agree_0_1.get(r.target, np.nan), repro_class=info.repro_class.get(r.target, "NA"), n_shared_true=int(n[r.target])))
for qcl, dcl in (("HepG2", "HepG2"), ("K562", "K562")): evaluate("same", qcl, dcl)
for qcl, dcl in (("HepG2", "K562"), ("K562", "HepG2")): evaluate("cross", qcl, dcl)
D = pd.DataFrame(out); D.to_csv(f"{O}/h3_encode_retrieval.tsv", sep="\t", index=False, float_format="%.4f")
S = []
for (st, qc), g in D.groupby(["setting", "query_cell"]):
    for cls, h in [("all", g)] + [(c, g[g.repro_class == c]) for c in ("high", "mid", "low")]:
        if len(h) == 0: continue
        N = h.n_candidates.median(); S.append(dict(setting=st, query_cell=qc, repro_class=cls, n_queries=len(h), n_rbps=h.target.nunique(), median_rank=h["rank"].median(), median_rank_frac=h.rank_frac.median(),
                      top1=(h["rank"] == 1).mean(), top5=(h["rank"] <= 5).mean(), top10=(h["rank"] <= 10).mean(), expected_top5=5 / N, expected_top10=10 / N, family_top5=h.family_top5.mean(), n_candidates=N))
S = pd.DataFrame(S); S.to_csv(f"{O}/h3_encode_summary.tsv", sep="\t", index=False, float_format="%.4f")
print("\n== H3 retrieval summary ==\n", S.round(3).to_string(index=False))
print("\n== most frequent wrong top-1 (attractors) ==")
for (st, qc), g in D.groupby(["setting", "query_cell"]): w = g[g["rank"] > 1].top1.value_counts().head(5); print(st, qc, w.to_dict())
print("\n== cross-cell, high-reproducibility RBPs ==\n", D[(D.setting == "cross") & (D.repro_class == "high")][["query_cell", "target", "assay", "n_query_events", "rank", "n_candidates", "score_true", "top1", "top5"]].sort_values("rank").round(3).to_string(index=False))
print("\n== same-cell leave-one-run-out ==\n", D[D.setting == "same"][["query_cell", "target", "assay", "n_query_events", "rank", "n_candidates", "score_true", "top1", "repro_class"]].sort_values("rank").round(3).to_string(index=False))
print("DONE", f"{O}/h3_encode_retrieval.tsv", f"{O}/h3_encode_summary.tsv")
