"""RBP x event matrices of RBP-mean control-centred dPSI on events covered in >= 80 % of RBPs in both cell lines.
Output: 13_model/programme/matrices.npz {events, targets_<cell>, Y_<cell> (float32, NaN = not measured), C_<cell> (bool, clean in >= 1 run)}."""
import numpy as np, pandas as pd, pyarrow.parquet as pq
PB = "/well/blagden/users/eic789/PerturbBridge"; cells = ("HepG2", "K562")
D = {c: pq.read_table(f"{PB}/14_inverse/dict_splicing_rbp_{c}.parquet").to_pandas().dropna(subset=["dpsi_cc"]) for c in cells}
cov = [D[c].groupby("event_id").size() / D[c].target.nunique() for c in cells]
ev = np.array(sorted(set(cov[0].index[cov[0] >= 0.8]) & set(cov[1].index[cov[1] >= 0.8]))); out = {"events": ev}
for c in cells:
    d = D[c][D[c].event_id.isin(ev)]; t = np.array(sorted(d.target.unique())); ti = pd.Index(t).get_indexer(d.target); ei = pd.Index(ev).get_indexer(d.event_id)
    Y = np.full((len(t), len(ev)), np.nan, np.float32); Y[ti, ei] = d.dpsi_cc.values; C = np.zeros(Y.shape, bool); C[ti, ei] = d.n_clean.values > 0
    out[f"targets_{c}"] = t; out[f"Y_{c}"] = Y; out[f"C_{c}"] = C
    print(c, Y.shape, "measured %.3f" % np.isfinite(Y).mean(), "clean pairs", int(C.sum()), flush=True)
np.savez_compressed(f"{PB}/13_model/programme/matrices.npz", **out); print("DONE")
