"""Predicted response-programme atlas for all human RNA-binding proteins (PLAN.md addendum b), and the post-freeze HeyA8 LARP1 knockout check.
Model: selected descriptor set (sequence/text + network), ridge, one model per cell line trained on all ENCODE RBPs of that cell line, reference
programme basis (as 04_annotate_programmes.py). Feature transforms are re-derived exactly as in 01_rbp_descriptors.py / 03_programme_model.py, checked
against the stored descriptors, and applied unchanged to proteins outside ENCODE.
Outputs: 13_model/programme/atlas_predicted_programmes.tsv (one row per protein and cell-line model), atlas_cv_<cell>.tsv (out-of-fold predictions for
ENCODE RBPs on the reference basis), inhouse_heya8_programme_check.tsv (in-house, not committed)."""
import numpy as np, pandas as pd, pickle, h5py, re, sys
PB = "/well/blagden/users/eic789/PerturbBridge"; O = f"{PB}/13_model/programme"; F = f"{PB}/06_features/rbp_functional"; A = f"{F}/atlas"; K = 10
def P(*x): print(*x, flush=True)
folds = pd.read_csv(f"{PB}/11_splits/rbp_folds.tsv", sep="\t"); T = folds.target.values; n = len(T); fold = folds.set_index("target").fold; D = np.load(f"{F}/rbp_descriptors.npz", allow_pickle=True); ids = pd.read_csv(f"{F}/rbp_id_map.tsv", sep="\t")
at = pd.read_csv(f"{A}/atlas_proteins.tsv", sep="\t"); new = at[~at.in_encode].reset_index(drop=True); m = len(new); P("ENCODE RBPs", n, "| new atlas proteins", m)
# ---- sequence/text block: transform fitted on the 478 ENCODE RBPs
R = {}
for c in ("HepG2", "K562"):
    z = np.load(f"{PB}/13_model/tensors_{c}.npz", allow_pickle=True)
    for t, a_, b_ in zip(z["targets"], z["R_esm"], z["R_gp"]): R[t] = np.hstack([a_, b_])
X0 = np.stack([R[t] for t in T]).astype(np.float64); mu0 = X0.mean(0); sd0 = X0.std(0) + 1e-6; X = (X0 - mu0) / sd0; mu1 = X.mean(0); Vt = np.linalg.svd(X - mu1, full_matrices=False)[2][:32]; seq = (X - mu1) @ Vt.T
P("seq transform reproduces stored descriptors: max |diff|", float(np.abs(np.abs(seq) - np.abs(D["seq"])).max()))
gp = pickle.load(open(f"{PB}/06_features/rbp_annotation/genept/GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle", "rb")); En = np.load(f"{A}/atlas_new_esm2.npy"); ei = pd.read_csv(f"{A}/atlas_new_esm2_index.tsv", sep="\t"); epos = {u: i for i, u in enumerate(ei.uniprot)}
has_gp = np.array([g in gp for g in new.gene]); Xn0 = np.stack([np.hstack([En[epos[u]], np.asarray(gp[g]) if g in gp else np.zeros(1536)]) for u, g in zip(new.uniprot, new.gene)]); seq_new = ((Xn0 - mu0) / sd0 - mu1) @ Vt.T
sgn = np.sign((seq * D["seq"]).sum(0)); seq_new = seq_new * sgn; seq = seq * sgn
# ---- network block: STRING embedding PCs over all proteins, Complex Portal membership projected on the ENCODE-fitted components
h = h5py.File(f"{F}/string/9606.protein.network.embeddings.v12.0.h5", "r"); E = h["embeddings"][:].astype(np.float32); pr = np.array([x.decode() if isinstance(x, bytes) else x for x in h["proteins"][:]]); ppos = {p: i for i, p in enumerate(pr)}
muE = E.mean(0); VE = np.linalg.svd(E - muE, full_matrices=False)[2][:32]; Ep = (E - muE) @ VE.T
cp = pd.read_csv(f"{F}/complexportal/9606.tsv", sep="\t"); upos = {u: i for i, u in enumerate(ids.uniprot)}; npos = {u: i for i, u in enumerate(new.uniprot)}; M = []; Mn = []
for parts in cp["Expanded participant list"].fillna(""):
    v = np.zeros(n, np.float32); w = np.zeros(m, np.float32)
    for u in {re.split(r"[-(]", p)[0] for p in parts.split("|")}:
        if u in upos: v[upos[u]] = 1
        if u in npos: w[npos[u]] = 1
    if v.sum() >= 2: M.append(v); Mn.append(w)
M = np.stack(M, 1); Mn = np.stack(Mn, 1); Mu, first = np.unique(M, axis=1, return_index=True); Mn = Mn[:, first]; muM = Mu.mean(0); VM = np.linalg.svd(Mu - muM, full_matrices=False)[2][:16]
net = np.zeros((n, 48)); has_net = np.zeros(n, bool)
for i, s in enumerate(ids.string_id):
    if s in ppos: net[i, :32] = Ep[ppos[s]]; has_net[i] = True
net[:, 32:] = (Mu - muM) @ VM.T; sg = np.sign((net * D["net"]).sum(0)); sg[sg == 0] = 1; net = net * sg
P("net transform reproduces stored descriptors: max |diff|", float(np.abs(net - D["net"]).max()))
net_new = np.zeros((m, 48)); has_net_new = np.zeros(m, bool)
for i, s in enumerate(new.string_id.fillna("")):
    if s in ppos: net_new[i, :32] = Ep[ppos[s]]; has_net_new[i] = True
net_new[:, 32:] = (Mn - muM) @ VM.T; net_new = net_new * sg; new["n_complexes"] = Mn.sum(1).astype(int)
P("new proteins: with STRING embedding", int(has_net_new.sum()), "| with GenePT text", int(has_gp.sum()), "| in >= 1 shared complex", int((Mn.sum(1) > 0).sum()))
blocks = {"seq": (seq, np.ones(n, bool), seq_new, np.ones(m, bool)), "net": (net, has_net, net_new, has_net_new)}
def design(idx):
    Xt, Xa = [], []
    for b in ("seq", "net"):
        Xb, hb, Xnb, hnb = blocks[b]; x = Xb[idx].copy(); h_ = hb[idx]; mu = x[h_].mean(0); sd = x[h_].std(0) + 1e-9; x = (x - mu) / sd; x[~h_] = 0; xa = np.vstack([(Xb - mu) / sd, (Xnb - mu) / sd]); ha = np.concatenate([hb, hnb]); xa[~ha] = 0
        if not h_.all(): x = np.hstack([x, ((h_ - h_.mean()) / (h_.std() + 1e-9))[:, None]]); xa = np.hstack([xa, ((ha - h_.mean()) / (h_.std() + 1e-9))[:, None]])
        Xt.append(x / np.sqrt(x.shape[1])); Xa.append(xa / np.sqrt(x.shape[1]))
    return np.hstack(Xt), np.hstack(Xa)
def ridge_fit(X, S, lam): mx = X.mean(0); ms = S.mean(0); Xc = X - mx; W = np.linalg.solve(Xc.T @ Xc + lam * np.eye(X.shape[1]), Xc.T @ (S - ms)); return lambda Z: (Z - mx) @ W + ms
LAMS = [0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 300]; z = np.load(f"{O}/matrices.npz", allow_pickle=True); ev = z["events"]; tpos = {t: i for i, t in enumerate(T)}; rows = []; basis = {}
names = np.concatenate([T, new.gene.values]); unip = np.concatenate([ids.uniprot.values, new.uniprot.values]); in_enc = np.concatenate([np.ones(n, bool), np.zeros(m, bool)]); hn_all = np.concatenate([has_net, has_net_new])
for c in ("HepG2", "K562"):
    Tc = z[f"targets_{c}"]; Y = z[f"Y_{c}"].astype(np.float64); mu = np.nanmean(Y, 0); Yc = np.where(np.isfinite(Y), Y, mu) - mu; V = np.linalg.svd(Yc, full_matrices=False)[2][:K]; S = Yc @ V.T; i0 = list(Tc).index("SF3B1")
    for k in range(K):
        if (k == 0 and S[i0, 0] < 0) or (k > 0 and np.abs(S[:, k]).max() != S[:, k].max()): V[k] *= -1; S[:, k] *= -1
    idx = np.array([tpos[t] for t in Tc]); Xt, Xa = design(idx); sdS = S.std(0); muS = S.mean(0); Sz = S / sdS; f = fold.reindex(Tc).values; err = []; oof = {}
    for lam in LAMS:
        pr_ = np.zeros_like(Sz)
        for k in range(5): te = f == k; pr_[te] = ridge_fit(Xt[~te], Sz[~te], lam)(Xt[te])
        err.append(((Sz - pr_) ** 2).sum()); oof[lam] = pr_
    lam = LAMS[int(np.argmin(err))]; cv = oof[lam] * sdS; pred = ridge_fit(Xt, Sz, lam)(Xa) * sdS; zp = (pred - muS) / sdS
    P(f"{c}: lambda {lam}; out-of-fold on the reference basis: programme-1 r {np.corrcoef(cv[:, 0], S[:, 0])[0, 1]:.3f}, spearman {pd.Series(cv[:, 0]).corr(pd.Series(S[:, 0]), method='spearman'):.3f}")
    pd.DataFrame(np.hstack([S, cv]), index=Tc, columns=[f"measured_P{k + 1}" for k in range(K)] + [f"cv_pred_P{k + 1}" for k in range(K)]).to_csv(f"{O}/atlas_cv_{c}.tsv", sep="\t", float_format="%.4f")
    Xn_ = Xa / (np.linalg.norm(Xa, axis=1, keepdims=True) + 1e-12); Xtn = Xt / (np.linalg.norm(Xt, axis=1, keepdims=True) + 1e-12); sim = Xn_ @ Xtn.T
    for i in range(len(names)):
        s_ = sim[i].copy()
        if in_enc[i] and names[i] in tpos and (idx == tpos[names[i]]).any(): s_[np.where(idx == tpos[names[i]])[0][0]] = -np.inf
        nb = np.argsort(-s_)[:3]; rows.append(dict(gene=names[i], uniprot=unip[i], model_cell_line=c, in_encode=bool(in_enc[i]), measured_in_this_cell_line=bool(in_enc[i] and names[i] in set(Tc)), has_string=bool(hn_all[i]), nearest_encode_rbps=",".join(Tc[nb]), nearest_similarity=round(float(s_[nb[0]]), 3), **{f"P{k + 1}_z": round(float(zp[i, k]), 3) for k in range(K)}))
    basis[c] = (mu, V, muS, sdS, zp)
At = pd.DataFrame(rows); At["P1_percentile_among_all"] = At.groupby("model_cell_line").P1_z.rank(pct=True).round(3); At.to_csv(f"{O}/atlas_predicted_programmes.tsv", sep="\t", index=False)
P("atlas rows", len(At), "| proteins", At.gene.nunique()); pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
for c in ("HepG2", "K562"):
    a_ = At[(At.model_cell_line == c) & ~At.in_encode].sort_values("P1_z", ascending=False); P(f"\n== {c} model: proteins outside ENCODE with the highest predicted programme-1 (core-spliceosome axis) score"); P(a_[["gene", "P1_z", "P2_z", "P3_z", "nearest_encode_rbps", "has_string"]].head(25).to_string(index=False))
P("\n== LARP family and selected proteins"); P(At[At.gene.isin(["LARP1", "LARP1B", "LARP4", "LARP4B", "LARP6", "LARP7", "SSB", "HNRNPA1", "PCBP2", "SF3B1"])].sort_values(["gene", "model_cell_line"])[["gene", "model_cell_line", "in_encode", "P1_z", "P2_z", "P3_z", "P4_z", "P5_z", "P1_percentile_among_all", "nearest_encode_rbps"]].to_string(index=False))
# ---- HeyA8 LARP1 knockout: project the measured knockout-minus-wild-type dPSI on the programmes (in-house, post-freeze)
sys.path.insert(0, f"{PB}/14_inverse"); import retrieval as Rt
Dh = "/well/blagden/users/eic789/LncRNASPLICE/output_lnc/star_salmon/rmats/KO-WT/rmats_post"; evp = pd.read_parquet(f"{Rt.U}/events.parquet", columns=["key", "event_id"]); k2i = dict(zip(evp.key, evp.event_id)); hr = []
for t, cols in Rt.COORD.items():
    x = pd.read_csv(f"{Dh}/{t}.MATS.JC.txt", sep="\t", dtype=str); key = t + "|" + x["chr"] + "|" + x["strand"] + "|" + x[cols].astype(str).agg("|".join, axis=1); x["event_id"] = key.map(k2i); x = x[x.event_id.notna()]
    c1 = x.IJC_SAMPLE_1.str.split(",", expand=True).astype(float).values + x.SJC_SAMPLE_1.str.split(",", expand=True).astype(float).values; c2 = x.IJC_SAMPLE_2.str.split(",", expand=True).astype(float).values + x.SJC_SAMPLE_2.str.split(",", expand=True).astype(float).values
    hr.append(pd.DataFrame({"event_id": x.event_id.values, "dpsi": x.IncLevelDifference.astype(float).values, "cov1": c1.mean(1), "cov2": c2.mean(1)}))
H = pd.concat(hr); H = H[(H.cov1 >= 10) & (H.cov2 >= 10) & np.isfinite(H.dpsi)].drop_duplicates("event_id").set_index("event_id").dpsi.reindex(ev); P(f"\n== HeyA8 LARP1 knockout: {int(H.notna().sum())} of {len(ev)} programme events measured")
chk = []
for c in ("HepG2", "K562"):
    mu, V, muS, sdS, zp = basis[c]; y = np.where(np.isfinite(H.values), H.values, mu) - mu; zm = (y @ V.T - muS) / sdS; cos = (zp @ zm) / (np.linalg.norm(zp, axis=1) * np.linalg.norm(zm) + 1e-12); il = int(np.where(names == "LARP1")[0][0]); rk = int((cos > cos[il]).sum()) + 1
    P(f" {c} model: measured HeyA8 z on P1..P10: " + " ".join(f"{v:+.2f}" for v in zm)); P(f"          predicted LARP1 z:              " + " ".join(f"{v:+.2f}" for v in zp[il]))
    P(f"          cosine(predicted LARP1, measured) {cos[il]:+.3f}; rank of LARP1 among {len(names)} atlas proteins {rk} (top {100 * rk / len(names):.1f} %); programme-1 sign match {bool(np.sign(zm[0]) == np.sign(zp[il, 0]))}; best-matching proteins: " + ", ".join(f"{names[j]}({cos[j]:.2f})" for j in np.argsort(-cos)[:8]))
    chk.append(dict(model_cell_line=c, cosine=cos[il], rank=rk, n=len(names), p1_measured_z=zm[0], p1_predicted_z=zp[il, 0], **{f"measured_P{k + 1}_z": zm[k] for k in range(K)}, **{f"predicted_P{k + 1}_z": zp[il, k] for k in range(K)}))
pd.DataFrame(chk).to_csv(f"{O}/inhouse_heya8_programme_check.tsv", sep="\t", index=False, float_format="%.4f"); P("DONE")
