"""RBP-level folds and per-run feature tiers for the forward model (design 2.4, 4.1).
Tier per run: A = the target has an IDR eCLIP in the same cell line; B = no same-cell eCLIP (C = ESM2/text only is an ablation
applied to every run, not a data tier). Folds: 5 folds by RBP, stratified by whether the RBP has any tier-A run, all runs of an
RBP (both assays, both cell lines) in one fold, and paralogs (BLASTP identity >= 40% over >= 30% of the shorter protein,
connected components, 11_splits/blast/allvsall.tsv) kept in the same fold. Family hold-out: RNA-binding-domain family from Pfam (first match in priority order).
Outputs: 11_splits/rbp_folds.tsv, run_tiers.tsv, leakage_report.txt, paralog_groups.tsv"""
import numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; O = f"{PB}/11_splits"; rng = np.random.default_rng(20261107)
pairs = pd.read_csv(f"{PB}/09_bam_layer/kd_control_pairs_all.tsv", sep="\t"); pairs["run"] = pairs.kd + "_vs_" + pairs.control
man = pd.read_csv(f"{PB}/01_encode_metadata/download_manifest.tsv", sep="\t", dtype=str)
ec = man[(man.assay == "eCLIP") & (man.bio_reps == "1;2")].groupby("target").cell_line.apply(set).to_dict()
pairs["tier"] = ["A" if r.cell_line in ec.get(r.target, set()) else "B" for r in pairs.itertuples()]
pairs["eclip_other_cell"] = [bool(ec.get(r.target, set()) - {r.cell_line}) for r in pairs.itertuples()]
pf = pd.read_csv(f"{PB}/06_features/rbp_annotation/rbp_uniprot_pfam.tsv", sep="\t")
FAM = [("PF00076", "RRM"), ("PF13893", "RRM"), ("PF14259", "RRM"), ("PF00013", "KH"), ("PF07650", "KH"), ("PF00642", "zf-CCCH"), ("PF00270", "DEAD"),
       ("PF04851", "DEAD"), ("PF00271", "DEAD"), ("PF00098", "zf-CCHC"), ("PF00575", "S1"), ("PF01423", "LSM"), ("PF00806", "PUF"), ("PF00096", "zf-C2H2"),
       ("PF00641", "zf-RanBP"), ("PF00035", "dsRBD"), ("PF00567", "TUDOR"), ("PF01424", "R3H"), ("PF05172", "Nup"), ("PF08064", "PWI"), ("PF05891", "Methyltransf")]
def fam(s):
    s = "" if pd.isna(s) else s
    for pid, name in FAM:
        if pid in s: return name
    return "other"
pf["family"] = pf.pfam.map(fam)
rbp = pairs.groupby("target").agg(n_runs=("run", "size"), cell_lines=("cell_line", lambda x: ",".join(sorted(set(x)))), assays=("assay", lambda x: ",".join(sorted(set(a.split()[0] for a in x)))),
                                 any_tier_A=("tier", lambda x: "A" in set(x))).reset_index().merge(pf[["target", "family", "pfam"]], on="target", how="left")
# paralog groups from all-vs-all BLASTP
bl = pd.read_csv(f"{O}/blast/allvsall.tsv", sep="\t", names=["q", "s", "pident", "length", "qlen", "slen", "evalue", "bits"])
bl = bl[(bl.q != bl.s) & (bl.pident >= 40) & (bl.length >= 0.3 * np.minimum(bl.qlen, bl.slen))].sort_values("bits", ascending=False).drop_duplicates(["q", "s"])
parent = {t: t for t in rbp.target}
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for r in bl.itertuples():
    if r.q in parent and r.s in parent: parent[find(r.q)] = find(r.s)
rbp["paralog_group"] = [find(t) for t in rbp.target]
grp = rbp.groupby("paralog_group").agg(n=("target", "size"), any_A=("any_tier_A", "any"), members=("target", lambda x: ",".join(sorted(x)))).reset_index()
grp["fold"] = -1
for a, g in grp.groupby("any_A"):
    idx = rng.permutation(g.index); grp.loc[idx, "fold"] = np.arange(len(idx)) % 5
grp.sort_values("n", ascending=False).to_csv(f"{O}/paralog_groups.tsv", sep="\t", index=False)
rbp = rbp.merge(grp[["paralog_group", "fold"]], on="paralog_group")
pairs = pairs.merge(rbp[["target", "fold", "family"]], on="target")
rbp.to_csv(f"{O}/rbp_folds.tsv", sep="\t", index=False); pairs.to_csv(f"{O}/run_tiers.tsv", sep="\t", index=False)
# leakage checks
fold_of = dict(zip(rbp.target, rbp.fold))
bla = pd.read_csv(f"{O}/blast/allvsall.tsv", sep="\t", names=["q", "s", "pident", "length", "qlen", "slen", "evalue", "bits"])
bla = bla[(bla.q != bla.s)].sort_values("bits", ascending=False).drop_duplicates(["q", "s"])
bla["cov"] = (bla.length / np.minimum(bla.qlen, bla.slen)).round(2); bla["q_fold"] = bla.q.map(fold_of); bla["s_fold"] = bla.s.map(fold_of)
cross = bla[bla.q_fold != bla.s_fold].sort_values("pident", ascending=False)
with open(f"{O}/leakage_report.txt", "w") as f:
    p = lambda *a: print(*a, file=f)
    p("RBPs:", len(rbp), "| runs:", len(pairs)); p("runs per tier:", pairs.tier.value_counts().to_dict(), "| RBPs with any tier-A run:", int(rbp.any_tier_A.sum()))
    p("tier A by cell line:", pairs[pairs.tier == "A"].groupby("cell_line").target.nunique().to_dict())
    p("fold sizes (RBPs):", rbp.fold.value_counts().sort_index().to_dict()); p("fold sizes (runs):", pairs.fold.value_counts().sort_index().to_dict())
    p("tier-A RBPs per fold:", rbp[rbp.any_tier_A].fold.value_counts().sort_index().to_dict())
    p("check: every RBP in exactly one fold:", pairs.groupby("target").fold.nunique().eq(1).all())
    p("check: every run in exactly one fold:", pairs.run.is_unique)
    ctrl = pairs.groupby("control").fold.nunique(); p("control experiments reused across folds:", int((ctrl > 1).sum()), "/", len(ctrl), "(control BAMs are unlabeled baselines; shared controls are not label leakage but are noted)")
    p("family counts:", rbp.family.value_counts().to_dict())
    p("families by fold:\n" + pd.crosstab(rbp.family, rbp.fold).to_string())
    p("\nparalog groups (BLASTP identity >= 40%, coverage >= 30%):", len(grp), "| groups with >1 member:", int((grp.n > 1).sum()), "| largest:", grp.sort_values("n", ascending=False).head(5)[["members"]].members.tolist())
    p("\nhighest-identity BLAST hits across folds (should all be < 40% identity or < 30% coverage):\n" + cross.head(15)[["q", "s", "pident", "cov", "evalue", "q_fold", "s_fold"]].to_string(index=False))
    p("cross-fold pairs with identity >= 40% and coverage >= 0.3:", int(((cross.pident >= 40) & (cross["cov"] >= 0.3)).sum()))
cross.to_csv(f"{O}/blast_cross_fold_hits.tsv", sep="\t", index=False)
print(open(f"{O}/leakage_report.txt").read())
