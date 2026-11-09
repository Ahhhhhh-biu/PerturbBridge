"""Gene-expression layer for the forward model: pack per cell line everything needed to predict the knockdown response of
gene expression (control-centred delta log2 TPM) for held-out RBPs, with the SAME runs / folds / tiers / RBP features as the
splicing tensors (13_model/tensors_<cell>.npz), so the two output layers are directly comparable.
  genes  : RSEM gene TPM (03_suppa_psi/matrices/<cell>_gene_TPM.tsv), mean TPM >= 1 over all samples; log2(TPM+1)
  target : Y[run, gene] = mean log2TPM of the KD replicates - leave-one-out consensus of the OTHER KD runs sharing the same
           control (the expression analogue of dpsi_cc; control-independent). Runs whose control is shared with no other KD are dropped.
  sig    : replicate-consistent flag: both KD replicates deviate from the consensus in the same direction, each by > 0.2,
           and |Y| > 0.5 (log2 units) -> used for loss weighting, model selection and sign metrics
  is_target : the knocked-down gene itself (gene_name == run target), excluded from training rows and metrics
  gene features Xg : log2 mean control TPM (all control samples of the cell line), log10 gene length, protein_coding flag,
           GenePT text embedding of the gene symbol (1536, zero + missing flag)
  eCLIP  : per (RBP with eCLIP in this cell line, gene): log1p(number of IDR peaks in the gene body, same strand), max signalValue
  RBP features : copied from tensors_<cell>.npz (ESM2, GenePT, motif), same target order
Also measurement references for this layer (13_model/expr_within_run.tsv, expr_ceiling_cross_run.tsv, expr_zero_point.tsv):
  within-run split half r(kd_rep1 - ctrl_rep1, kd_rep2 - ctrl_rep2); cross-run same-RBP r of Y and sign agreement on sig genes of run A
  scored in run B; zero point = random pairs of different RBPs. Needs bedtools on PATH.
Usage: python 10_build_expr_tensors.py <cell_line> -> 13_model/tensors_expr_<cell_line>.npz"""
import sys, gzip, pickle, subprocess, numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; cl = sys.argv[1]; M = f"{PB}/13_model"; TMP = f"{PB}/envs/tmp"; rng = np.random.default_rng(0)
tiers = pd.read_csv(f"{PB}/11_splits/run_tiers.tsv", sep="\t"); tiers = tiers[tiers.cell_line == cl].set_index("run")
z = np.load(f"{M}/tensors_{cl}.npz", allow_pickle=True); targets = z["targets"]; tpos = {x: i for i, x in enumerate(targets)}
assert set(tiers.target) <= set(targets), "run targets not all in splicing tensors"
# genes
G = pd.read_csv(f"{PB}/03_suppa_psi/matrices/{cl}_gene_TPM.tsv", sep="\t", index_col=0); G = G[G.mean(axis=1) >= 1.0]; L = np.log2(G + 1)
ann = []
for line in gzip.open(f"{PB}/03_suppa_psi/gencode.v29.primary_assembly.annotation.gtf.gz", "rt"):
    if line.startswith("#"): continue
    f = line.split("\t", 9)
    if f[2] != "gene": continue
    a = dict(x.strip().split(" ", 1) for x in f[8].strip().strip(";").split(";") if x.strip()); ann.append((a["gene_id"].strip('"'), a["gene_name"].strip('"'), a["gene_type"].strip('"'), f[0], int(f[3]), int(f[4]), f[6]))
ann = pd.DataFrame(ann, columns=["gene_id", "gene_name", "gene_type", "chr", "start", "end", "strand"]).set_index("gene_id")
L = L.loc[L.index.isin(ann.index)]; genes = np.array(L.index); ann = ann.loc[genes]; gpos = {g: i for i, g in enumerate(genes)}
def reps(acc): return [c for c in L.columns if c.startswith(acc + "_rep")]
runs = [r for r in tiers.index if len(reps(tiers.kd[r])) >= 2 and len(reps(tiers.control[r])) >= 1]
kd_mean = pd.DataFrame({r: L[reps(tiers.kd[r])].mean(axis=1) for r in runs}); ctrl_mean = pd.DataFrame({r: L[reps(tiers.control[r])].mean(axis=1) for r in runs})
kd1 = pd.DataFrame({r: L[reps(tiers.kd[r])[0]] for r in runs}); kd2 = pd.DataFrame({r: L[reps(tiers.kd[r])[1]] for r in runs})
ctrl = tiers.control[runs]; cons = {}
for r in runs:
    others = [o for o in runs if ctrl[o] == ctrl[r] and o != r]
    if others: cons[r] = kd_mean[others].mean(axis=1)
runs = [r for r in runs if r in cons]; cons = pd.DataFrame(cons)[runs]
Y = (kd_mean[runs] - cons).T; D1 = (kd1[runs] - cons).T; D2 = (kd2[runs] - cons).T; RAW = (kd_mean[runs] - ctrl_mean[runs]).T   # runs x genes
SIG = (np.sign(D1.values) == np.sign(D2.values)) & (np.minimum(np.abs(D1.values), np.abs(D2.values)) > 0.2) & (np.abs(Y.values) > 0.5)
IS_T = np.array([[ann.gene_name[g] == tiers.target[r] for g in genes] for r in runs])
print(cl, "genes", len(genes), "runs", len(runs), "of", len(tiers), "| sig frac", SIG.mean().round(4), "per run median", int(np.median(SIG.sum(1))), "| Y sd", Y.values.std().round(3), "| target gene found in", IS_T.any(1).sum(), "runs")
# gene features
ctrl_cols = [c for c in L.columns if c.rsplit("_rep", 1)[0] in set(tiers.control)]
gp = pickle.load(open(f"{PB}/06_features/rbp_annotation/genept/GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle", "rb"))
G_gp = np.array([np.asarray(gp[n]) if n in gp else np.zeros(1536) for n in ann.gene_name], dtype=np.float32); G_has = np.array([n in gp for n in ann.gene_name], dtype=np.float32)
Xg = np.stack([L[ctrl_cols].mean(axis=1).values, np.log10(ann.end - ann.start + 1).values, (ann.gene_type == "protein_coding").values.astype(float), G_has], 1).astype(np.float32)
Xg_names = np.array(["log2_ctrl_tpm", "log10_len", "protein_coding", "has_genept"])
# eCLIP peaks in gene bodies
man = pd.read_csv(f"{PB}/01_encode_metadata/download_manifest.tsv", sep="\t", dtype=str); pk = man[(man.assay == "eCLIP") & (man.bio_reps == "1;2") & (man.cell_line == cl) & man.target.isin(tpos)]
with open(f"{TMP}/_expr_peaks_{cl}.bed", "w") as out:
    for r in pk.itertuples():
        for line in gzip.open(f"{PB}/02_encode_download/{r.dest}", "rt"):
            f = line.rstrip("\n").split("\t"); out.write(f"{f[0]}\t{f[1]}\t{f[2]}\t{r.target}\t0\t{f[5]}\t{f[6]}\n")
pd.DataFrame({"chr": ann.chr, "s": ann.start - 1, "e": ann.end, "name": genes, "sc": 0, "strand": ann.strand}).to_csv(f"{TMP}/_expr_genes_{cl}.bed", sep="\t", header=False, index=False)
cmd = f"bedtools intersect -a {TMP}/_expr_genes_{cl}.bed -b {TMP}/_expr_peaks_{cl}.bed -s -wa -wb | cut -f4,10,13 > {TMP}/_expr_hits_{cl}.tsv"
subprocess.run(["bash", "-c", "set -o pipefail; " + cmd], check=True)
h = pd.read_csv(f"{TMP}/_expr_hits_{cl}.tsv", sep="\t", header=None, names=["gene", "rbp", "sig"]); ec_targets = np.array(sorted(h.rbp.unique())); ecpos = {x: i for i, x in enumerate(ec_targets)}
agg = h.groupby(["rbp", "gene"]).sig.agg(["size", "max"]).reset_index(); EC = np.zeros((len(ec_targets), len(genes), 2), dtype=np.float32)
EC[[ecpos[x] for x in agg.rbp], [gpos[g] for g in agg.gene]] = np.stack([np.log1p(agg["size"].values), agg["max"].values], 1)
target_ecidx = np.array([ecpos.get(x, -1) for x in targets]); subprocess.run(["bash", "-c", f"rm -f {TMP}/_expr_*_{cl}*"])
print("eCLIP RBPs with peaks in gene bodies:", len(ec_targets), "| genes with any peak:", int((EC[..., 0] > 0).any(0).sum()), "| mean genes bound per RBP:", int((EC[..., 0] > 0).sum(1).mean()))
ri = tiers.loc[runs]; np.savez_compressed(f"{M}/tensors_expr_{cl}.npz", genes=genes, gene_name=ann.gene_name.values, Xg=Xg, Xg_names=Xg_names, G_gp=G_gp, targets=targets,
    R_esm=z["R_esm"], R_gp=z["R_gp"], R_mot=z["R_mot"], R_has_mot=z["R_has_mot"], ec_targets=ec_targets, EC=EC, target_ecidx=target_ecidx,
    runs=np.array(runs), run_target=np.array([tpos[x] for x in ri.target]), run_assay=ri.assay.str.startswith("CRISPR").values.astype(np.float32), run_tier=(ri.tier == "A").values.astype(np.float32),
    run_fold=ri.fold.values, run_control=ri.control.values, Y=Y.values.astype(np.float32), RAW=RAW.values.astype(np.float32), SIG=SIG, IS_T=IS_T)
# measurement references
W = []
for r in runs:
    c = reps(tiers.control[r]); c1, c2 = (L[c[0]], L[c[1]]) if len(c) >= 2 else (L[c[0]], L[c[0]])
    d1 = kd1[r] - c1; d2 = kd2[r] - c2; W.append(dict(cell_line=cl, run=r, target=tiers.target[r], assay=tiers.assay[r], n_sig=int(SIG[runs.index(r)].sum()), r_split_half=np.corrcoef(d1, d2)[0, 1], ctrl_reps=len(c)))
W = pd.DataFrame(W); W.to_csv(f"{M}/expr_within_run_{cl}.tsv", sep="\t", index=False)
X, Z = [], []; rt = tiers.target[runs]; by = pd.Series(range(len(runs)), index=runs)
for tg, rr in rt.groupby(rt).groups.items():
    rr = list(rr)
    if len(rr) < 2: continue
    for ra in rr:
        for rb in rr:
            if ra == rb: continue
            ia, ib = by[ra], by[rb]; s = SIG[ia] & ~IS_T[ia] & ~IS_T[ib]
            X.append(dict(cell_line=cl, target=tg, run_select=ra, assay_select=tiers.assay[ra], run_test=rb, assay_test=tiers.assay[rb], n_sig_a=int(s.sum()), r_all=np.corrcoef(Y.values[ia], Y.values[ib])[0, 1],
                          r_sig=np.corrcoef(Y.values[ia][s], Y.values[ib][s])[0, 1] if s.sum() > 2 else np.nan, sign_sig=(np.sign(Y.values[ia][s]) == np.sign(Y.values[ib][s])).mean() if s.sum() else np.nan))
for _ in range(400):
    ia, ib = rng.choice(len(runs), 2, replace=False)
    if rt.iloc[ia] == rt.iloc[ib]: continue
    s = SIG[ia] & ~IS_T[ia]; Z.append(dict(cell_line=cl, run_a=runs[ia], run_b=runs[ib], r_all=np.corrcoef(Y.values[ia], Y.values[ib])[0, 1], sign_sig=(np.sign(Y.values[ia][s]) == np.sign(Y.values[ib][s])).mean() if s.sum() else np.nan))
X = pd.DataFrame(X); Z = pd.DataFrame(Z); X.to_csv(f"{M}/expr_ceiling_cross_run_{cl}.tsv", sep="\t", index=False); Z.to_csv(f"{M}/expr_zero_point_{cl}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(f"\n== {cl} within-run split-half r (kd_rep - ctrl_rep): mean {W.r_split_half.mean():.3f} median {W.r_split_half.median():.3f} | by RBP >0.3: {(W.groupby('target').r_split_half.mean() > 0.3).sum()} / {W.target.nunique()}")
print(f"== {cl} cross-run same RBP ({len(X)} ordered pairs): r_all mean {X.r_all.mean():.3f} | r_sig {X.r_sig.mean():.3f} | sign_sig {X.sign_sig.mean():.3f} | cross-technology r_all {X[X.assay_select != X.assay_test].r_all.mean():.3f}")
print(X.round(3).to_string(index=False))
print(f"== {cl} zero point, different RBPs ({len(Z)} pairs): r_all mean {Z.r_all.mean():.3f} | sign_sig {Z.sign_sig.mean():.3f}")
print("DONE", f"{M}/tensors_expr_{cl}.npz")
