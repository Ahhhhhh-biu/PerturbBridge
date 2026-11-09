"""Control for reagent-specific (same shRNA hairpin) off-target signal in the cross-cell retrieval: the dictionary of the other cell line
is built from runs of ONE technology only, and queries are split by their own technology. Cells: shRNA->shRNA (same hairpin possible),
CRISPR->CRISPR (same guide possible), shRNA->CRISPR and CRISPR->shRNA (no shared reagent). If retrieval survives across technologies the
signal is biological. Output 14_inverse/h3_cross_technology.tsv (per query) and summary to stdout (logs/h3_cross_technology.txt)."""
import sys, numpy as np, pandas as pd; sys.path.insert(0, "/well/blagden/users/eic789/PerturbBridge/14_inverse"); import retrieval as R
O = R.O; pd.set_option("display.width", 250); info = R.rbp_info(); out = []
runs = {cl: pd.read_parquet(f"{O}/dict_splicing_runs_{cl}.parquet") for cl in ("HepG2", "K562")}
prof = {}
for cl, t in runs.items():
    for asy in ("shRNA RNA-seq", "CRISPR RNA-seq"):
        g = t[t.assay == asy].groupby(["event_id", "target"]).dpsi_cc.mean().reset_index(); prof[(cl, asy)] = g.pivot(index="event_id", columns="target", values="dpsi_cc"); print(cl, asy, prof[(cl, asy)].shape)
for qcl, dcl in (("HepG2", "K562"), ("K562", "HepG2")):
    t = runs[qcl]; tq = t[t.clean]; rt = t[["run", "target", "assay"]].drop_duplicates()
    for r in rt.itertuples():
        q = tq[tq.run == r.run].set_index("event_id").dpsi_cc
        if len(q) < 50: continue
        for dasy in ("shRNA RNA-seq", "CRISPR RNA-seq"):
            P = prof[(dcl, dasy)]
            if r.target not in P.columns: continue
            s, n = R.score_against(q, P); s = s.dropna()
            if r.target not in s.index: continue
            rank = int((s > s[r.target]).sum() + 1); top = s.sort_values(ascending=False)
            out.append(dict(query_cell=qcl, query_assay=r.assay.split()[0], dict_assay=dasy.split()[0], same_technology=r.assay == dasy, run=r.run, target=r.target, n_query=len(q), n_candidates=len(s), rank=rank, score_true=s[r.target], top1=top.index[0], repro_class=info.repro_class.get(r.target, "NA")))
D = pd.DataFrame(out); D.to_csv(f"{O}/h3_cross_technology.tsv", sep="\t", index=False, float_format="%.4f")
S = D.groupby(["query_assay", "dict_assay", "same_technology"]).agg(n_queries=("run", "size"), n_rbps=("target", "nunique"), n_candidates=("n_candidates", "median"), median_rank=("rank", "median"), top1=("rank", lambda x: (x == 1).mean()), top5=("rank", lambda x: (x <= 5).mean()), top10=("rank", lambda x: (x <= 10).mean())).reset_index()
print("\n== cross-cell retrieval by technology pairing (all RBPs) ==\n", S.round(3).to_string(index=False))
S2 = D.groupby(["same_technology", "repro_class"]).agg(n_queries=("run", "size"), median_rank=("rank", "median"), top1=("rank", lambda x: (x == 1).mean()), top5=("rank", lambda x: (x <= 5).mean()), top10=("rank", lambda x: (x <= 10).mean())).reset_index()
print("\n== by reproducibility class ==\n", S2.round(3).to_string(index=False))
# paired within query run: same-technology dictionary rank vs cross-technology dictionary rank for the same query
pv = D.pivot_table(index="run", columns="same_technology", values="rank"); pv = pv.dropna(); print(f"\nqueries with both dictionaries: {len(pv)} | median rank same-tech {pv[True].median():.0f} vs cross-tech {pv[False].median():.0f} | top5 same {np.mean(pv[True] <= 5):.3f} vs cross {np.mean(pv[False] <= 5):.3f} | cross-tech top5 for high-repro: {np.mean(D[(~D.same_technology) & (D.repro_class == 'high')]['rank'] <= 5):.3f} (n={int(((~D.same_technology) & (D.repro_class == 'high')).sum())})")
print("DONE")
