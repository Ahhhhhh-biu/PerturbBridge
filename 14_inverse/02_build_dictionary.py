"""Block 1: the RBP response dictionary used by the inverse mapping (H3/H4).
Measured ENCODE knockdown profiles only (model-predicted profiles are not used: the forward model carries transferable signal only
for the core spliceosome class, STATUS 2026-11-08/09). Per cell line:
  dict_splicing_runs_<cell>.parquet : one row per (run, event) with control-centred dPSI (dpsi_cc) and the replicate-consistent 'clean' flag
  dict_splicing_rbp_<cell>.parquet  : per (target, event): mean dpsi_cc over the target's runs, n_runs, n_clean
  dict_expression_rbp_<cell>.parquet: per (target, gene): mean control-centred delta log2 TPM, n_runs, n_sig
  dict_rbps.tsv                      : per (target, cell line): n_runs, within-run replicate direction agreement (reproducibility),
                                       class high (>0.7) / mid (0.6-0.7) / low, n clean events, expression split-half r
Usage: python 02_build_dictionary.py"""
import numpy as np, pandas as pd, pyarrow.parquet as pq
PB = "/well/blagden/users/eic789/PerturbBridge"; U = f"{PB}/10_event_universe/encode_full"; O = f"{PB}/14_inverse"
rep = pd.read_csv(f"{PB}/13_model/ceiling_within_run_by_rbp.tsv", sep="\t")[["cell_line", "target", "agree_0_1", "n_runs"]]
R = []
for cl in ("HepG2", "K562"):
    t = pq.read_table(f"{U}/train_table_{cl}.parquet", columns=["event_id", "run", "target", "assay", "dpsi_cc"]).to_pandas()
    t = t[np.isfinite(t.dpsi_cc)].reset_index(drop=True); z = np.load(f"{PB}/13_model/tensors_{cl}.npz", allow_pickle=True)
    assert len(t) == len(z["row_y"]) and np.allclose(t.dpsi_cc.values, z["row_y"], atol=1e-5), "row order differs from tensors"
    t["clean"] = np.load(f"{PB}/13_model/clean_flags_{cl}.npz")["row_clean"]; t["dpsi_cc"] = t.dpsi_cc.astype(np.float32)
    t[["event_id", "run", "target", "assay", "dpsi_cc", "clean"]].to_parquet(f"{O}/dict_splicing_runs_{cl}.parquet", index=False)
    g = t.groupby(["target", "event_id"], sort=False).agg(dpsi_cc=("dpsi_cc", "mean"), n_runs=("run", "size"), n_clean=("clean", "sum")).reset_index()
    g["dpsi_cc"] = g.dpsi_cc.astype(np.float32); g["n_runs"] = g.n_runs.astype(np.int8); g["n_clean"] = g.n_clean.astype(np.int8); g.to_parquet(f"{O}/dict_splicing_rbp_{cl}.parquet", index=False)
    nclean = t[t.clean].groupby("target").event_id.nunique().rename("n_clean_events"); nruns = t.groupby("target").run.nunique().rename("n_runs_dict")
    print(cl, "rows", len(t), "| targets", t.target.nunique(), "| rbp x event", len(g), "| clean events per target median", int(nclean.median()))
    # expression
    ze = np.load(f"{PB}/13_model/tensors_expr_{cl}.npz", allow_pickle=True); Y = ze["Y"]; SIG = ze["SIG"]; IS_T = ze["IS_T"]; tg = ze["targets"][ze["run_target"]]
    rows = []
    for x in np.unique(tg):
        m = tg == x; y = np.where(IS_T[m], np.nan, Y[m]); rows.append(pd.DataFrame({"target": x, "gene_id": ze["genes"], "gene_name": ze["gene_name"], "dlog2tpm_cc": np.nanmean(y, 0).astype(np.float32), "n_runs": np.int8(m.sum()), "n_sig": SIG[m].sum(0).astype(np.int8)}))
    ex = pd.concat(rows, ignore_index=True); ex = ex[np.isfinite(ex.dlog2tpm_cc)]; ex.to_parquet(f"{O}/dict_expression_rbp_{cl}.parquet", index=False)
    w = pd.read_csv(f"{PB}/13_model/expr_within_run_{cl}.tsv", sep="\t").groupby("target").r_split_half.mean().rename("expr_split_half_r")
    d = rep[rep.cell_line == cl].set_index("target").join(nclean).join(nruns).join(w); d["in_expression"] = d.index.isin(ex.target.unique())
    d["repro_class"] = np.select([d.agree_0_1 > 0.7, d.agree_0_1 > 0.6], ["high", "mid"], "low"); R.append(d.reset_index())
    print(cl, "expression rbp x gene", len(ex), "| classes", d.repro_class.value_counts().to_dict())
D = pd.concat(R); D.to_csv(f"{O}/dict_rbps.tsv", sep="\t", index=False, float_format="%.4f")
both = D.groupby("target").cell_line.nunique(); print("targets in both cell lines:", int((both == 2).sum()), "| high in either:", D[D.repro_class == "high"].target.nunique(), "| high in both:", int((D[D.repro_class == "high"].groupby("target").size() == 2).sum()))
print("DONE", O)
