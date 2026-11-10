"""RBP descriptor blocks for the programme-level model (all unsupervised with respect to ENCODE splicing data).
seq: ESM2 + GenePT (32 PCs over the 478 RBPs); ps: K562 genome-wide Perturb-seq pseudo-bulk signature (Replogle et al. 2022, normalised bulk; 32 PCs
fitted on all perturbations, own gene column zeroed, plus signature norm); net: STRING v12 network embedding (32 PCs over all proteins) + Complex Portal
membership (16 SVD components); dep: DepMap 24Q4 CRISPR gene effect profile across cell lines (32 PCs over all genes).
Output: 06_features/rbp_functional/rbp_descriptors.npz {targets, <block>, has_<block>} and rbp_id_map.tsv."""
import numpy as np, pandas as pd, h5py, gzip, re
PB = "/well/blagden/users/eic789/PerturbBridge"; F = f"{PB}/06_features/rbp_functional"
def pcs(X, k, fit=None):
    fit = X if fit is None else fit; mu = fit.mean(0); U, S, Vt = np.linalg.svd(fit - mu, full_matrices=False); return ((X - mu) @ Vt[:k].T).astype(np.float32)
folds = pd.read_csv(f"{PB}/11_splits/rbp_folds.tsv", sep="\t"); T = folds.target.values; n = len(T); out = {"targets": T}
# identifiers: target -> UniProt -> STRING protein -> preferred name, Ensembl gene, Entrez
up = pd.read_csv(f"{PB}/06_features/rbp_annotation/rbp_uniprot_pfam.tsv", sep="\t").set_index("target").uniprot.reindex(T)
al = pd.read_csv(f"{F}/string/9606.protein.aliases.v12.0.txt.gz", sep="\t"); al.columns = ["sid", "alias", "source"]
u2s = al[al.source.str.contains("UniProt_AC") & al.alias.isin(set(up.dropna()))].drop_duplicates("alias").set_index("alias").sid
sid = up.map(u2s); info = pd.read_csv(f"{F}/string/9606.protein.info.v12.0.txt.gz", sep="\t").set_index("#string_protein_id").preferred_name
ensg = al[al.source == "Ensembl_gene"].drop_duplicates("sid").set_index("sid").alias; entrez = al[al.source == "Ensembl_HGNC_entrez_id"].drop_duplicates("sid").set_index("sid").alias
ids = pd.DataFrame({"target": T, "uniprot": up.values, "string_id": sid.values, "symbol": sid.map(info).values, "ensg": sid.map(ensg).values, "entrez": sid.map(entrez).values})
# seq block from the training tensors
R = {}
for c in ("HepG2", "K562"):
    z = np.load(f"{PB}/13_model/tensors_{c}.npz", allow_pickle=True)
    for t, a, b in zip(z["targets"], z["R_esm"], z["R_gp"]): R[t] = np.hstack([a, b])
X = np.stack([R[t] for t in T]); X = (X - X.mean(0)) / (X.std(0) + 1e-6); out["seq"] = pcs(X, 32); out["has_seq"] = np.ones(n, bool)
# Perturb-seq
g = h5py.File(f"{F}/perturbseq/K562_gwps_normalized_bulk_01.h5ad", "r"); P = g["X"][:].astype(np.float32); P[~np.isfinite(P)] = 0
obs = np.array([x.decode() for x in g["obs"][g["obs"].attrs.get("_index", "index")][:]]); var = np.array([x.decode() for x in g["var"][g["var"].attrs.get("_index", "index")][:]])
pg = np.array([o.split("_")[-1] for o in obs]); psym = np.array([o.split("_")[1] for o in obs]); vpos = {v: i for i, v in enumerate(var)}
for i, e in enumerate(pg):
    if e in vpos: P[i, vpos[e]] = 0          # own-gene knockdown readout removed
norm = np.sqrt((P ** 2).mean(1)); Pp = pcs(P, 32); ps = np.zeros((n, 33), np.float32); has = np.zeros(n, bool)
for i, (t, e, s) in enumerate(zip(T, ids.ensg, ids.symbol)):
    m = (pg == e) | (psym == t) | (psym == s)
    if m.any(): ps[i, :32] = Pp[m].mean(0); ps[i, 32] = norm[m].mean(); has[i] = True
out["ps"] = ps; out["has_ps"] = has; ids["perturbseq_rows"] = [int(((pg == e) | (psym == t) | (psym == s)).sum()) for t, e, s in zip(T, ids.ensg, ids.symbol)]
# STRING network embedding + Complex Portal
h = h5py.File(f"{F}/string/9606.protein.network.embeddings.v12.0.h5", "r"); E = h["embeddings"][:].astype(np.float32); pr = np.array([x.decode() if isinstance(x, bytes) else x for x in h["proteins"][:]]); ppos = {p: i for i, p in enumerate(pr)}
Ep = pcs(E, 32); net = np.zeros((n, 48), np.float32); has = np.zeros(n, bool)
for i, s in enumerate(ids.string_id):
    if s in ppos: net[i, :32] = Ep[ppos[s]]; has[i] = True
cp = pd.read_csv(f"{F}/complexportal/9606.tsv", sep="\t"); upos = {u: i for i, u in enumerate(ids.uniprot)}; M = []
for parts in cp["Expanded participant list"].fillna(""):
    v = np.zeros(n, np.float32)
    for u in {re.split(r"[-(]", p)[0] for p in parts.split("|")}:
        if u in upos: v[upos[u]] = 1
    if v.sum() >= 2: M.append(v)
M = np.unique(np.stack(M, 1), axis=1); U, S, Vt = np.linalg.svd(M - M.mean(0), full_matrices=False); net[:, 32:] = (U[:, :16] * S[:16]).astype(np.float32)
out["net"] = net; out["has_net"] = has; ids["n_complexes"] = M.sum(1).astype(int); print("complexes with >= 2 RBPs", M.shape[1], "| RBPs in >= 1", int((M.sum(1) > 0).sum()))
# DepMap
d = pd.read_csv(f"{F}/depmap/CRISPRGeneEffect.csv", index_col=0); G = d.values.T.astype(np.float32); G = np.where(np.isfinite(G), G, np.nanmean(G, 1, keepdims=True)); G = G - G.mean(0)
gent = np.array([re.search(r"\((\d+)\)", c).group(1) for c in d.columns]); gsym = np.array([c.split(" ")[0] for c in d.columns]); Gp = pcs(G, 32); dep = np.zeros((n, 32), np.float32); has = np.zeros(n, bool)
for i, (t, e, s) in enumerate(zip(T, ids.entrez, ids.symbol)):
    m = (gent == str(e)) | (gsym == t) | (gsym == s)
    if m.any(): dep[i] = Gp[m][0]; has[i] = True
out["dep"] = dep; out["has_dep"] = has
for b in ("seq", "ps", "net", "dep"): ids[f"has_{b}"] = out[f"has_{b}"]; print(b, out[b].shape, "coverage", int(out[f"has_{b}"].sum()), "of", n)
ids.to_csv(f"{F}/rbp_id_map.tsv", sep="\t", index=False); np.savez_compressed(f"{F}/rbp_descriptors.npz", **out); print("DONE")
