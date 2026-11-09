#!/usr/bin/env python
"""H2 go/no-go: is the knockdown response more perturbation-specific at the splicing level (ΔPSI)
than at the expression level (Δlog2 TPM)?

Per KD experiment e with control c (same cell line, ENCODE possible_controls):
  Δ_e = mean(reps of e) − mean(reps of c)        (per event / per gene)
Three statistics, each computed for ΔPSI and Δlog2TPM on the same experiments:
  (1) perturbation-specific variance ratio per feature: Var_e(Δ_e) / mean_e(Var_within_e),
      where Var_within_e = Var(rep1 − rep2)/2 (replicate noise of the KD arm)
  (2) leave-one-out mean-baseline R²: how well the mean Δ of all other KDs predicts Δ_e
      (high = a mean baseline is hard to beat; Ahlmann-Eltze / Csendes critique)
  (3) replicate-split self-consistency: corr(Δ_e using rep1, Δ_e using rep2) — the measurement ceiling
Outputs: h2_summary.tsv, h2_per_experiment.tsv, h2_figure.png
"""
import pandas as pd, numpy as np, os, sys
PB = "/well/blagden/users/eic789/PerturbBridge"; OUT = f"{PB}/04_h2_check"
ex = pd.read_csv(f"{PB}/01_encode_metadata/experiments.tsv", sep="\t", dtype=str).fillna("")
ex = ex[ex.cell_line.isin(["K562", "HepG2"]) & ex.assay.isin(["shRNA RNA-seq", "CRISPR RNA-seq"])]
MIN_PSI_SAMPLES = 0.8   # event must have PSI defined in ≥80% of samples
MIN_TPM = 1.0           # gene mean TPM filter

def reps(M, acc):
    c = [x for x in M.columns if x.startswith(acc + "_rep")]
    return M[c]

rows, per = [], {}
for cl in ["K562", "HepG2"]:
    psi = pd.read_csv(f"{PB}/03_suppa_psi/psi/{cl}_events.psi", sep="\t", index_col=0)
    sheet = pd.read_csv(f"{PB}/03_suppa_psi/matrices/{cl}_transcript_samples.tsv", sep="\t", dtype=str)
    psi = psi[[c for c in psi.columns if c in set(sheet[sheet.quant_tool == "rsem"]["sample"])]]   # RSEM only (each replicate also has a kallisto file)
    psi = psi[psi.notna().mean(axis=1) >= MIN_PSI_SAMPLES]
    g = pd.read_csv(f"{PB}/03_suppa_psi/matrices/{cl}_gene_TPM.tsv", sep="\t", index_col=0)
    g = np.log2(g[g.mean(axis=1) >= MIN_TPM] + 1)
    for label, M in [("dPSI", psi), ("dlog2TPM", g)]:
        D, D1, D2, W = {}, {}, {}, {}
        for e in ex[(ex.cell_line == cl) & (ex.control_type == "")].itertuples():
            kd = reps(M, e.accession)
            ctrls = [reps(M, c) for c in e.controls.split(";") if c]
            if kd.shape[1] < 2 or not ctrls: continue
            ctrl = pd.concat(ctrls, axis=1)
            if ctrl.shape[1] < 1: continue
            d = kd.mean(axis=1) - ctrl.mean(axis=1)
            D[e.accession] = d
            D1[e.accession] = kd.iloc[:, 0] - ctrl.mean(axis=1); D2[e.accession] = kd.iloc[:, 1] - ctrl.mean(axis=1)
            W[e.accession] = ((kd.iloc[:, 0] - kd.iloc[:, 1]) ** 2) / 2
        D = pd.DataFrame(D); D1 = pd.DataFrame(D1); D2 = pd.DataFrame(D2); W = pd.DataFrame(W)
        ok = D.notna().mean(axis=1) >= 0.8
        D, D1, D2, W = D[ok], D1[ok], D2[ok], W[ok]
        # (1) perturbation-specific variance ratio per feature
        ratio = D.var(axis=1, skipna=True) / W.mean(axis=1, skipna=True).replace(0, np.nan)
        # (2) leave-one-out mean-baseline R² per experiment
        S = D.sum(axis=1, skipna=True); n = D.notna().sum(axis=1)
        loo_r2, selfc = {}, {}
        for a in D.columns:
            m = (S - D[a].fillna(0)) / (n - D[a].notna().astype(int)).replace(0, np.nan)
            v = D[a].notna() & m.notna()
            y, yhat = D.loc[v, a], m[v]
            loo_r2[a] = 1 - ((y - yhat) ** 2).sum() / ((y - y.mean()) ** 2).sum()
            vv = D1[a].notna() & D2[a].notna()
            selfc[a] = np.corrcoef(D1.loc[vv, a], D2.loc[vv, a])[0, 1]
        per[(cl, label)] = pd.DataFrame({"loo_mean_baseline_R2": loo_r2, "replicate_split_corr": selfc})
        rows.append({"cell_line": cl, "readout": label, "n_experiments": D.shape[1], "n_features": D.shape[0],
                     "median_pert_specific_var_ratio": float(np.nanmedian(ratio)),
                     "frac_features_ratio_gt_2": float((ratio > 2).mean()),
                     "median_loo_mean_baseline_R2": float(np.nanmedian(list(loo_r2.values()))),
                     "median_replicate_split_corr": float(np.nanmedian(list(selfc.values())))})
        D.to_csv(f"{OUT}/{cl}_{label}_delta_matrix.tsv.gz", sep="\t", float_format="%.4f")
        print(rows[-1], file=sys.stderr)
summary = pd.DataFrame(rows); summary.to_csv(f"{OUT}/h2_summary.tsv", sep="\t", index=False)
pe = pd.concat({k: v for k, v in per.items()}, names=["cell_line", "readout", "experiment"]).reset_index()
pe.to_csv(f"{OUT}/h2_per_experiment.tsv", sep="\t", index=False)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(7, 3), dpi=300)
for i, col in enumerate(["loo_mean_baseline_R2", "replicate_split_corr"]):
    data = [pe[(pe.cell_line == cl) & (pe.readout == r)][col].dropna() for cl in ["K562", "HepG2"] for r in ["dPSI", "dlog2TPM"]]
    ax[i].boxplot(data, tick_labels=[f"{cl}\n{r}" for cl in ["K562", "HepG2"] for r in ["ΔPSI", "Δlog2TPM"]], showfliers=False)
    ax[i].set_ylabel(col.replace("_", " ")); ax[i].tick_params(labelsize=7)
plt.tight_layout(); plt.savefig(f"{OUT}/h2_figure.png")
print(summary.to_string(), file=sys.stderr)
