#!/usr/bin/env python
"""Junction-level (rMATS JC) test of the splicing target on the 12 shRNA/CRISPR same-RBP pairs.
Per KD-vs-control run: ΔPSI = IncLevelDifference per event (SE, A3SS, A5SS, MXE, RI), coverage filter on
mean(IJC+SJC) per group. Reports (i) cross-technology r of ΔPSI for the same RBP, vs random-pair null,
stratified by coverage; (ii) replicate-split ceiling using per-replicate IncLevel (rep1 vs rep1, rep2 vs rep2);
(iii) number of events with FDR<0.05 & |ΔPSI|>0.1 per run."""
import pandas as pd, numpy as np, glob, os, itertools, sys
PB = "/well/blagden/users/eic789/PerturbBridge"; D = f"{PB}/08_junction_psi"
pairs = pd.read_csv(f"{D}/kd_control_pairs.tsv", sep="\t")
TYPES = ["SE", "A3SS", "A5SS", "MXE", "RI"]
def load_run(kd, ct, min_cov=10):
    out = []
    for t in TYPES:
        p = f"{D}/rmats/{kd}_vs_{ct}/{t}.MATS.JC.txt"
        if not os.path.exists(p): continue
        d = pd.read_csv(p, sep="\t", low_memory=False)
        key_cols = [c for c in ["chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
                                "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE",
                                "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "riExonStart_0base", "riExonEnd"] if c in d.columns]
        d["event"] = t + ":" + d[key_cols].astype(str).agg(":".join, axis=1)
        def covm(col): return d[col].astype(str).str.split(",").apply(lambda v: np.mean([float(x) for x in v if x not in ("NA", "nan")]) if any(x not in ("NA", "nan") for x in v) else np.nan)
        d["cov1"] = covm("IJC_SAMPLE_1") + covm("SJC_SAMPLE_1"); d["cov2"] = covm("IJC_SAMPLE_2") + covm("SJC_SAMPLE_2")
        def repl(col, i): return d[col].astype(str).str.split(",").apply(lambda v: float(v[i]) if len(v) > i and v[i] not in ("NA", "nan") else np.nan)
        d["dpsi_r1"] = repl("IncLevel1", 0) - repl("IncLevel2", 0); d["dpsi_r2"] = repl("IncLevel1", 1) - repl("IncLevel2", 1)
        out.append(d[["event", "IncLevelDifference", "FDR", "cov1", "cov2", "dpsi_r1", "dpsi_r2"]])
    r = pd.concat(out).set_index("event")
    return r[(r.cov1 >= min_cov) & (r.cov2 >= min_cov)]
def corr(a, b):
    v = a.notna() & b.notna(); return (np.corrcoef(a[v], b[v])[0, 1], int(v.sum())) if v.sum() > 30 else (np.nan, int(v.sum()))
runs = {}
for r in pairs.itertuples():
    if os.path.exists(f"{D}/rmats/{r.kd}_vs_{r.control}/SE.MATS.JC.txt"): runs[r.kd] = load_run(r.kd, r.control)
print("runs loaded:", len(runs), file=sys.stderr)
rows_sig, rows_x = [], []
for r in pairs.itertuples():
    if r.kd not in runs: continue
    d = runs[r.kd]; c, n = corr(d.dpsi_r1, d.dpsi_r2)
    rows_sig.append({"kd": r.kd, "target": r.target, "assay": r.assay, "cell_line": r.cell_line, "n_events_cov10": len(d),
                     "n_sig_fdr05_dpsi0.1": int(((d.FDR < 0.05) & (d.IncLevelDifference.abs() > 0.1)).sum()), "replicate_split_r": c})
meta = pairs.set_index("kd"); rng = np.random.default_rng(0)
for cl in pairs.cell_line.unique():
    sub = pairs[pairs.cell_line == cl]
    for t in sub.target.unique():
        sh = sub[(sub.target == t) & (sub.assay == "shRNA RNA-seq")].kd.tolist(); cr = sub[(sub.target == t) & (sub.assay == "CRISPR RNA-seq")].kd.tolist()
        for a, b in itertools.product(sh, cr):
            if a not in runs or b not in runs: continue
            A, B = runs[a], runs[b]; common = A.index.intersection(B.index)
            c_all, n_all = corr(A.loc[common, "IncLevelDifference"], B.loc[common, "IncLevelDifference"])
            hi = common[(A.loc[common, "cov1"] >= 50) & (B.loc[common, "cov1"] >= 50)]
            c_hi, n_hi = corr(A.loc[hi, "IncLevelDifference"], B.loc[hi, "IncLevelDifference"])
            sig = common[((A.loc[common, "FDR"] < 0.05) & (A.loc[common, "IncLevelDifference"].abs() > 0.1)) | ((B.loc[common, "FDR"] < 0.05) & (B.loc[common, "IncLevelDifference"].abs() > 0.1))]
            c_sig, n_sig = corr(A.loc[sig, "IncLevelDifference"], B.loc[sig, "IncLevelDifference"])
            others = [k for k in runs if meta.loc[k, "cell_line"] == cl and meta.loc[k, "target"] != t and meta.loc[k, "assay"] == "CRISPR RNA-seq"]
            null = []
            for o in rng.choice(others, size=min(20, len(others)), replace=True):
                O = runs[o]; cm = A.index.intersection(O.index); null.append(corr(A.loc[cm, "IncLevelDifference"], O.loc[cm, "IncLevelDifference"])[0])
            rows_x.append({"cell_line": cl, "target": t, "shRNA": a, "CRISPR": b, "r_all": c_all, "n_all": n_all, "r_cov50": c_hi, "n_cov50": n_hi,
                           "r_sig_either": c_sig, "n_sig_either": n_sig, "null_median_r": float(np.nanmedian(null)) if null else np.nan})
S = pd.DataFrame(rows_sig); X = pd.DataFrame(rows_x)
S.to_csv(f"{D}/junction_per_run.tsv", sep="\t", index=False); X.to_csv(f"{D}/junction_cross_technology.tsv", sep="\t", index=False)
print("\n== per run ==", file=sys.stderr); print(S.round(3).to_string(), file=sys.stderr)
print("\n== shRNA vs CRISPR same RBP, junction-level ΔPSI ==", file=sys.stderr); print(X.round(3).to_string(), file=sys.stderr)
