"""Describe the response programmes: reference basis per cell line (SVD of the centred RBP x event matrix over all RBPs), variance explained,
RBPs with the largest scores, event-type composition of the top-loading events, relation of programme-1 score to the size of the knockdown
effect, and correspondence of programmes between the two cell lines. Output: 13_model/programme/programme_annotation.txt, programme_scores_<cell>.tsv"""
import numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; O = f"{PB}/13_model/programme"; z = np.load(f"{O}/matrices.npz", allow_pickle=True); ev = z["events"]; K = 10
rep = pd.read_csv(f"{PB}/14_inverse/dict_rbps.tsv", sep="\t").set_index(["cell_line", "target"]); out = open(f"{O}/programme_annotation.txt", "w"); V = {}
def P(*x): print(*x, flush=True); print(*x, file=out, flush=True)
et = np.array([e.split("_")[0] for e in ev])
for c in ("HepG2", "K562"):
    T = z[f"targets_{c}"]; Y = z[f"Y_{c}"].astype(np.float64); mu = np.nanmean(Y, 0); Yc = np.where(np.isfinite(Y), Y, mu) - mu; U, S, Vt = np.linalg.svd(Yc, full_matrices=False)
    i = list(T).index("SF3B1") if "SF3B1" in T else 0; sc = Yc @ Vt[:K].T
    for k in range(K):
        if (k == 0 and sc[i, 0] < 0) or (k > 0 and np.abs(sc[:, k]).max() != sc[:, k].max()): Vt[k] *= -1; sc[:, k] *= -1
    V[c] = Vt[:K]; zs = (sc - sc.mean(0)) / sc.std(0); pd.DataFrame(sc, index=T, columns=[f"P{k + 1}" for k in range(K)]).to_csv(f"{O}/programme_scores_{c}.tsv", sep="\t", float_format="%.4f")
    nclean = rep.loc[c].n_clean_events.reindex(T).values; cls = rep.loc[c].repro_class.reindex(T).values; eff = np.sqrt((Yc ** 2).mean(1))
    P(f"\n== {c}: {len(T)} RBPs x {len(ev)} events; variance explained by programmes 1-10: " + " ".join(f"{v:.3f}" for v in (S[:K] ** 2 / (S ** 2).sum())) + f" (sum {(S[:K] ** 2).sum() / (S ** 2).sum():.3f})")
    P(f"   programme-1 score vs RMS effect size r {np.corrcoef(sc[:, 0], eff)[0, 1]:.2f}; vs number of clean events r {np.corrcoef(sc[:, 0], nclean)[0, 1]:.2f}; mean z by reproducibility class: " + ", ".join(f"{g} {zs[cls == g, 0].mean():+.2f}" for g in ("high", "mid", "low")))
    for k in range(K):
        o = np.argsort(-zs[:, k]); top = np.argsort(-np.abs(Vt[k]))[:500]; comp = pd.Series(et[top]).value_counts(normalize=True); base = pd.Series(et).value_counts(normalize=True)
        P(f" P{k + 1}: + " + ", ".join(f"{T[j]}({zs[j, k]:.1f})" for j in o[:8]) + " | - " + ", ".join(f"{T[j]}({zs[j, k]:.1f})" for j in o[::-1][:4]))
        P(f"      top-500 events: " + ", ".join(f"{t} {comp.get(t, 0):.2f} (all {base[t]:.2f})" for t in base.index) + f"; fraction with positive loading {np.mean(Vt[k][top] > 0):.2f}; RBPs with |z| > 3: {int((np.abs(zs[:, k]) > 3).sum())}")
Cc = np.abs(V["HepG2"] @ V["K562"].T); P("\n== correspondence between cell lines: |cosine| of loadings, rows HepG2 P1-P10, columns K562 P1-P10"); P(pd.DataFrame(Cc, index=[f"H{k + 1}" for k in range(K)], columns=[f"K{k + 1}" for k in range(K)]).round(2).to_string()); P("DONE")
