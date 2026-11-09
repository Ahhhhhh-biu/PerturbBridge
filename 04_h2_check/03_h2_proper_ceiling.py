#!/usr/bin/env python
"""Replicate-split ceiling WITHOUT shared control noise: D1 = kd_rep1 − ctrl_rep1, D2 = kd_rep2 − ctrl_rep2
(first listed control, 2 reps). ΔPSI and Δlog2TPM, overall and by event expression bin; reproducible responders recounted."""
import pandas as pd, numpy as np, sys
PB = "/well/blagden/users/eic789/PerturbBridge"; OUT = f"{PB}/04_h2_check"
ex = pd.read_csv(f"{PB}/01_encode_metadata/experiments.tsv", sep="\t", dtype=str).fillna("")
ex = ex[ex.cell_line.isin(["K562", "HepG2"]) & ex.assay.isin(["shRNA RNA-seq", "CRISPR RNA-seq"])]
ioe = pd.read_csv(f"{PB}/03_suppa_psi/events/gencode_v29_all_strict.ioe", sep="\t").set_index("event_id")
def reps(M, acc): return M[[c for c in M.columns if c.startswith(acc + "_rep")]]
def corr(a, b):
    v = a.notna() & b.notna(); return np.corrcoef(a[v], b[v])[0, 1] if v.sum() > 30 else np.nan
rows = []
for cl in ["K562", "HepG2"]:
    psi = pd.read_csv(f"{PB}/03_suppa_psi/psi/{cl}_events.psi", sep="\t", index_col=0)
    sheet = pd.read_csv(f"{PB}/03_suppa_psi/matrices/{cl}_transcript_samples.tsv", sep="\t", dtype=str)
    psi = psi[[c for c in psi.columns if c in set(sheet[sheet.quant_tool == "rsem"]["sample"])]]
    psi = psi[psi.notna().mean(axis=1) >= 0.8]
    g = pd.read_csv(f"{PB}/03_suppa_psi/matrices/{cl}_gene_TPM.tsv", sep="\t", index_col=0); g = np.log2(g[g.mean(axis=1) >= 1.0] + 1)
    T = pd.read_csv(f"{PB}/03_suppa_psi/matrices/{cl}_transcript_TPM_rsem_suppa.tsv", sep="\t", header=0, index_col=0)
    ctrl_cols = [c for c in T.columns if c.split("_rep")[0] in set(ex[(ex.cell_line == cl) & (ex.control_type != "")].accession)]
    tmean = T[ctrl_cols].mean(axis=1)
    ev_expr = pd.Series({e: tmean.reindex(str(ioe.loc[e, "total_transcripts"]).split(",")).sum() for e in psi.index})
    expr_bin = pd.cut(ev_expr, [-1, 1, 5, 20, 1e9], labels=["<1", "1-5", "5-20", ">20"]).astype(str)
    meta = ex.set_index("accession")
    for e in ex[(ex.cell_line == cl) & (ex.control_type == "")].itertuples():
        c0 = e.controls.split(";")[0]
        for label, M in [("dPSI", psi), ("dlog2TPM", g)]:
            kd, ct = reps(M, e.accession), reps(M, c0)
            if kd.shape[1] < 2 or ct.shape[1] < 2: continue
            D1, D2 = kd.iloc[:, 0] - ct.iloc[:, 0], kd.iloc[:, 1] - ct.iloc[:, 1]
            r = {"cell_line": cl, "experiment": e.accession, "assay": e.assay, "target": e.target, "readout": label,
                 "ceiling_r_all": corr(D1, D2),
                 "n_reproducible_responders": int(((D1.abs() > (0.1 if label == "dPSI" else 1.0)) & (D2.abs() > (0.1 if label == "dPSI" else 1.0)) & (np.sign(D1) == np.sign(D2))).sum())}
            if label == "dPSI":
                for lvl in ["<1", "1-5", "5-20", ">20"]:
                    idx = expr_bin.index[expr_bin == lvl]; r[f"ceiling_r_expr_{lvl}"] = corr(D1.loc[idx], D2.loc[idx])
            rows.append(r)
    print(cl, "done", file=sys.stderr)
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/h2c_proper_ceiling.tsv", sep="\t", index=False)
print(R.groupby(["cell_line", "readout"])[[c for c in R.columns if c.startswith("ceiling_r") or c.startswith("n_rep")]].median().round(3).to_string(), file=sys.stderr)
