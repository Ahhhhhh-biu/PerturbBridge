"""One 5-mer weight vector per RBP for the motif layer of the perturbation encoder (design 2.4 tiers A/B).
Source priority: (1) ENCODE RBNS 5-mer enrichment R (06_features/rbns/rbns_5mer_R.tsv): log2(R) clipped at 0, z-scored over the
1024 k-mers; (2) ATtRACT human, non-mutated motifs (06_features/attract/raw/ATtRACT_db.txt): IUPAC consensus expanded to all
5-mers contained in (motif length >= 5) or containing (length 4) the motif, weight = count of motif entries supporting the
5-mer, z-scored. Motif density of event region g for RBP r = <kmer_counts(g), vec(r)> / n_kmers(g), computed downstream from
10_event_universe/encode_full/kmer5_counts.npz.
Output: 06_features/rbp_annotation/rbp_motif_5mer.tsv (target, source, 1024 columns) and rbp_motif_sources.tsv"""
import itertools, numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; kmers = ["".join(p) for p in itertools.product("ACGT", repeat=5)]; kpos = {k: i for i, k in enumerate(kmers)}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
def z(v): v = np.asarray(v, float); return (v - v.mean()) / (v.std() + 1e-9)
vec, src = {}, []
R = pd.read_csv(f"{PB}/06_features/rbns/rbns_5mer_R.tsv", sep="\t", index_col=0)
for tg, row in R.iterrows():
    v = np.zeros(1024); 
    for k, x in row.items(): v[kpos[k.rstrip("*")]] = max(0.0, np.log2(max(x, 1e-6)))  # some tables mark significant k-mers with "*"
    vec[tg] = z(v); src.append((tg, "RBNS", int((row > 2).sum()), row.idxmax().rstrip("*")))
at = pd.read_csv(f"{PB}/06_features/attract/raw/ATtRACT_db.txt", sep="\t"); at = at[(at.Organism == "Homo_sapiens") & (at.Mutated == "no")]
for tg, g in at.groupby("Gene_name"):
    if tg in vec: continue
    v = np.zeros(1024); n = 0
    for m in g.Motif.astype(str).str.upper().unique():
        if any(c not in IUPAC for c in m): continue
        exp = ["".join(p) for p in itertools.product(*[IUPAC[c] for c in m])]
        hits = set()
        for s in exp:
            if len(s) >= 5: hits.update(s[i:i + 5] for i in range(len(s) - 4))
            elif len(s) == 4: hits.update(k for k in kmers if s in k)
        for k in hits: v[kpos[k]] += 1
        n += 1
    if v.sum() == 0: continue
    vec[tg] = z(v); src.append((tg, "ATtRACT", n, kmers[int(np.argmax(v))]))
M = pd.DataFrame(vec, index=kmers).T; M.index.name = "target"; S = pd.DataFrame(src, columns=["target", "source", "n_motifs_or_enriched", "top_5mer"])
M.insert(0, "source", S.set_index("target").loc[M.index, "source"].values); M.to_csv(f"{PB}/06_features/rbp_annotation/rbp_motif_5mer.tsv", sep="\t"); S.to_csv(f"{PB}/06_features/rbp_annotation/rbp_motif_sources.tsv", sep="\t", index=False)
kd = set(pd.read_csv(f"{PB}/09_bam_layer/kd_control_pairs_all.tsv", sep="\t").target); ec = set(pd.read_csv(f"{PB}/01_encode_metadata/download_manifest.tsv", sep="\t", dtype=str).query("assay=='eCLIP'").target)
print("RBPs with motif vector:", len(M), S.source.value_counts().to_dict(), "| in KD set:", len(set(M.index) & kd), "| KD without eCLIP but with motif:", len((set(M.index) & kd) - ec))
print(S[S.target.isin(["CELF1", "RBFOX2", "PTBP1", "HNRNPC", "SRSF1", "QKI", "TIA1", "U2AF2", "PUM1", "KHSRP"])].to_string(index=False))
