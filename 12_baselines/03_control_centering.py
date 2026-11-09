"""Shared-control artefact: dPSI = PSI(KD) - PSI(control); every KD compared against the same control experiment inherits that
control's replicate noise, so runs sharing a control correlate (null r ~0.49 vs ~0.02 for different controls, 01_ceiling.py).
Control-centred dPSI: for each run, subtract the leave-one-out mean dPSI of the other runs that use the same control
(per event; runs whose control is used by no other run keep the raw dPSI — nothing to centre against, and no sharing means no leakage). Writes <U>/train_table_<cl>.parquet with a new column dpsi_cc and reports the
null correlations and the cross-assay pairs before/after.
Usage: python 03_control_centering.py"""
import numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; U = f"{PB}/10_event_universe/encode_full"; rng = np.random.default_rng(2)
tiers = pd.read_csv(f"{PB}/11_splits/run_tiers.tsv", sep="\t").set_index("run")
def r(a, b):
    m = np.isfinite(a) & np.isfinite(b); return np.corrcoef(a[m], b[m])[0, 1] if m.sum() > 2 else np.nan
for cl in ["HepG2", "K562"]:
    t = pd.read_parquet(f"{U}/train_table_{cl}.parquet"); M = t.pivot(index="event_id", columns="run", values="dpsi")
    info = tiers.loc[M.columns]; ctrl = info.control.values; CC = M.copy()
    nper = pd.Series(ctrl).value_counts(); print(cl, "runs per control: n controls", len(nper), "| singletons", int((nper == 1).sum()), "| median", nper.median(), "| max", nper.max())
    for c in nper.index:
        cols = np.where(ctrl == c)[0]
        if len(cols) == 1: continue
        sub = M.iloc[:, cols].values; s = np.nansum(sub, 1, keepdims=True); n = np.isfinite(sub).sum(1, keepdims=True)
        loo = (s - np.where(np.isfinite(sub), sub, 0)) / np.maximum(n - 1, 1); loo = np.where((n - 1) == 0, np.nan, loo)
        CC.iloc[:, cols] = sub - loo
    cc = CC.stack().rename("dpsi_cc").reset_index(); t = t.drop(columns=[c for c in t.columns if c == "dpsi_cc"]).merge(cc, on=["event_id", "run"], how="left")
    t.to_parquet(f"{U}/train_table_{cl}.parquet", index=False)
    runs = list(M.columns); rows = []
    for _ in range(3000):
        a, b = rng.choice(len(runs), 2, replace=False)
        if info.target.values[a] == info.target.values[b]: continue
        rows.append((ctrl[a] == ctrl[b], r(M.values[:, a], M.values[:, b]), r(CC.values[:, a], CC.values[:, b])))
    nl = pd.DataFrame(rows, columns=["same_control", "r_raw", "r_cc"]); print(nl.groupby("same_control")[["r_raw", "r_cc"]].describe(percentiles=[.5, .9]).round(3).T.to_string())
    for tg, g in info.groupby("target"):
        if g.assay.nunique() == 2:
            a, b = g.index[0], g.index[1]; print(f"cross-assay {tg}: raw r={r(M[a].values, M[b].values):.3f}  centred r={r(CC[a].values, CC[b].values):.3f}")
    print("control-centred dPSI available for", int(CC.notna().values.sum()), "of", int(M.notna().values.sum()), "covered pairs")
