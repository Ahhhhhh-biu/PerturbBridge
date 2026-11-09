#!/usr/bin/env python
"""H2 follow-up: where is the reproducible ΔPSI signal?
 (a) replicate-split ceiling by event type and by event expression (control TPM of the event's transcripts)
 (b) ceiling restricted to responsive events (|ΔPSI| > 0.1) and number of reproducible responders per KD, by assay
 (c) shRNA vs CRISPR, same RBP, same cell line: cross-technology ΔPSI / Δlog2TPM correlation vs random-pair null
Outputs: h2b_*.tsv and h2b_figure.png"""
import pandas as pd, numpy as np, sys, itertools
PB = "/well/blagden/users/eic789/PerturbBridge"; OUT = f"{PB}/04_h2_check"
ex = pd.read_csv(f"{PB}/01_encode_metadata/experiments.tsv", sep="\t", dtype=str).fillna("")
ex = ex[ex.cell_line.isin(["K562", "HepG2"]) & ex.assay.isin(["shRNA RNA-seq", "CRISPR RNA-seq"])]
ioe = pd.read_csv(f"{PB}/03_suppa_psi/events/gencode_v29_all_strict.ioe", sep="\t").set_index("event_id")
def reps(M, acc): return M[[c for c in M.columns if c.startswith(acc + "_rep")]]
def corr(a, b):
    v = a.notna() & b.notna()
    return np.corrcoef(a[v], b[v])[0, 1] if v.sum() > 30 else np.nan
rows_a, rows_b, rows_c = [], [], []
for cl in ["K562", "HepG2"]:
    psi = pd.read_csv(f"{PB}/03_suppa_psi/psi/{cl}_events.psi", sep="\t", index_col=0)
    sheet = pd.read_csv(f"{PB}/03_suppa_psi/matrices/{cl}_transcript_samples.tsv", sep="\t", dtype=str)
    psi = psi[[c for c in psi.columns if c in set(sheet[sheet.quant_tool == "rsem"]["sample"])]]
    psi = psi[psi.notna().mean(axis=1) >= 0.8]
    g = pd.read_csv(f"{PB}/03_suppa_psi/matrices/{cl}_gene_TPM.tsv", sep="\t", index_col=0)
    g = np.log2(g[g.mean(axis=1) >= 1.0] + 1)
    # event expression: mean control TPM of the event's total transcripts
    T = pd.read_csv(f"{PB}/03_suppa_psi/matrices/{cl}_transcript_TPM_rsem_suppa.tsv", sep="\t", header=0, index_col=0)
    ctrl_cols = [c for c in T.columns if c.split("_rep")[0] in set(ex[(ex.cell_line == cl) & (ex.control_type != "")].accession)]
    tmean = T[ctrl_cols].mean(axis=1)
    ev_expr = pd.Series({e: tmean.reindex(str(ioe.loc[e, "total_transcripts"]).split(",")).sum() for e in psi.index})
    ev_type = pd.Series({e: e.split(";")[1].split(":")[0] for e in psi.index})
    expr_bin = pd.cut(ev_expr, [-1, 1, 5, 20, 1e9], labels=["<1", "1-5", "5-20", ">20"])
    D, D1, D2, G = {}, {}, {}, {}
    for e in ex[(ex.cell_line == cl) & (ex.control_type == "")].itertuples():
        kd = reps(psi, e.accession); ctrls = [reps(psi, c) for c in e.controls.split(";") if c]
        if kd.shape[1] < 2 or not ctrls: continue
        cm = pd.concat(ctrls, axis=1).mean(axis=1)
        D[e.accession] = kd.mean(axis=1) - cm; D1[e.accession] = kd.iloc[:, 0] - cm; D2[e.accession] = kd.iloc[:, 1] - cm
        kg = reps(g, e.accession); cg = pd.concat([reps(g, c) for c in e.controls.split(";") if c], axis=1).mean(axis=1)
        G[e.accession] = kg.mean(axis=1) - cg
    D, D1, D2, G = map(pd.DataFrame, (D, D1, D2, G))
    meta = ex.set_index("accession")
    # (a) ceiling by type / expression bin
    for a in D.columns:
        for name, grp in [("type", ev_type), ("expr", expr_bin.astype(str))]:
            for lvl in grp.unique():
                idx = grp.index[grp == lvl]
                rows_a.append({"cell_line": cl, "experiment": a, "strat": name, "level": lvl, "n_events": len(idx),
                               "replicate_split_corr": corr(D1.loc[idx, a], D2.loc[idx, a])})
    # (b) responsive events
    for a in D.columns:
        resp = D.index[(D[a].abs() > 0.1)]
        rep_resp = int(((D1[a].abs() > 0.1) & (D2[a].abs() > 0.1) & (np.sign(D1[a]) == np.sign(D2[a]))).sum())
        rows_b.append({"cell_line": cl, "experiment": a, "assay": meta.loc[a, "assay"], "target": meta.loc[a, "target"],
                       "n_resp_events_abs_gt_0.1": len(resp), "n_reproducible_responders": rep_resp,
                       "ceiling_corr_all": corr(D1[a], D2[a]), "ceiling_corr_responsive": corr(D1.loc[resp, a], D2.loc[resp, a])})
    # (c) shRNA vs CRISPR same RBP
    sh = meta[(meta.cell_line == cl) & (meta.assay == "shRNA RNA-seq") & (meta.control_type == "")]
    cr = meta[(meta.cell_line == cl) & (meta.assay == "CRISPR RNA-seq") & (meta.control_type == "")]
    common = sorted(set(sh.target) & set(cr.target))
    rng = np.random.default_rng(0)
    for t in common:
        for a, b in itertools.product(sh.index[sh.target == t], cr.index[cr.target == t]):
            if a not in D or b not in D: continue
            null_psi = [corr(D[a], D[rng.choice(cr.index[cr.target != t])]) for _ in range(20)]
            null_g = [corr(G[a], G[rng.choice(cr.index[cr.target != t])]) for _ in range(20)]
            rows_c.append({"cell_line": cl, "target": t, "shRNA": a, "CRISPR": b,
                           "dPSI_cross_tech_corr": corr(D[a], D[b]), "dPSI_null_median": float(np.nanmedian(null_psi)),
                           "dlog2TPM_cross_tech_corr": corr(G[a], G[b]), "dlog2TPM_null_median": float(np.nanmedian(null_g))})
    print(cl, "done", file=sys.stderr)
A = pd.DataFrame(rows_a); B = pd.DataFrame(rows_b); C = pd.DataFrame(rows_c)
A.to_csv(f"{OUT}/h2b_ceiling_stratified.tsv", sep="\t", index=False); B.to_csv(f"{OUT}/h2b_responsive_per_experiment.tsv", sep="\t", index=False); C.to_csv(f"{OUT}/h2b_cross_technology.tsv", sep="\t", index=False)
print("\n== (a) median replicate-split ceiling by stratum ==", file=sys.stderr)
print(A.groupby(["cell_line", "strat", "level"]).replicate_split_corr.median().round(3).to_string(), file=sys.stderr)
print("\n== (b) per-experiment medians by assay ==", file=sys.stderr)
print(B.groupby(["cell_line", "assay"])[["n_resp_events_abs_gt_0.1", "n_reproducible_responders", "ceiling_corr_all", "ceiling_corr_responsive"]].median().round(3).to_string(), file=sys.stderr)
print("\n== (c) shRNA vs CRISPR same RBP ==", file=sys.stderr); print(C.round(3).to_string(), file=sys.stderr)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(10, 3), dpi=300)
for lvl in ["<1", "1-5", "5-20", ">20"]:
    ax[0].plot([lvl], [A[(A.strat == "expr") & (A.level == lvl)].replicate_split_corr.median()], "o")
ax[0].set_title("ΔPSI ceiling by event expression (ctrl TPM)", fontsize=7); ax[0].set_ylabel("median replicate-split r", fontsize=7)
ax[1].scatter(B.n_reproducible_responders, B.ceiling_corr_responsive, s=4); ax[1].set_xscale("symlog"); ax[1].set_xlabel("reproducible responders per KD", fontsize=7); ax[1].set_ylabel("ceiling r (responsive events)", fontsize=7)
if len(C): ax[2].scatter(C.dPSI_null_median, C.dPSI_cross_tech_corr, s=8, label="ΔPSI"); ax[2].scatter(C.dlog2TPM_null_median, C.dlog2TPM_cross_tech_corr, s=8, label="Δlog2TPM"); ax[2].plot([-.2, 1], [-.2, 1], "k--", lw=.5); ax[2].set_xlabel("random-pair null r", fontsize=7); ax[2].set_ylabel("shRNA vs CRISPR same-RBP r", fontsize=7); ax[2].legend(fontsize=6)
for a in ax: a.tick_params(labelsize=6)
plt.tight_layout(); plt.savefig(f"{OUT}/h2b_figure.png")
