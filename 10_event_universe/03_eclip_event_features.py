"""Event-level eCLIP features on the coordinate-keyed event universe (replaces 06_features/01_eclip_event_features.py,
which was keyed on SUPPA event ids). For every event three regions are intersected with each RBP's IDR-merged eCLIP
narrowPeak (per cell line), same strand only:
  seg = the alternative segment (seg_start/seg_end from event_features.parquet)
  up  = 300 nt immediately upstream of the segment in transcript orientation
  dn  = 300 nt immediately downstream of the segment in transcript orientation
Output: <universe_dir>/eclip_event_features.parquet, one row per (event_id, rbp, cell_line) with at least one peak:
  sig_up, sig_seg, sig_dn = max narrowPeak signalValue (log2 fold enrichment) in the region, 0 if no peak
  nlp_up, nlp_seg, nlp_dn = max -log10 p in the region, 0 if no peak
Needs bedtools on PATH (module load BEDTools/2.31.0-GCC-12.3.0).
Usage: python 03_eclip_event_features.py <universe_dir>"""
import sys, gzip, subprocess, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; U = sys.argv[1]; W = 300
ev = pd.read_parquet(f"{U}/events.parquet", columns=["event_id", "chr", "strand"]).merge(
    pd.read_parquet(f"{U}/event_features.parquet", columns=["event_id", "seg_start", "seg_end"]), on="event_id")
plus = ev.strand == "+"; lo = (ev.seg_start - W).clip(lower=0); hi = ev.seg_end + W
reg = pd.concat([
    pd.DataFrame({"chr": ev.chr, "s": ev.seg_start, "e": ev.seg_end, "name": ev.event_id + "|seg", "sc": 0, "strand": ev.strand}),
    pd.DataFrame({"chr": ev.chr, "s": lo, "e": ev.seg_start, "name": ev.event_id + plus.map({True: "|up", False: "|dn"}), "sc": 0, "strand": ev.strand}),
    pd.DataFrame({"chr": ev.chr, "s": ev.seg_end, "e": hi, "name": ev.event_id + plus.map({True: "|dn", False: "|up"}), "sc": 0, "strand": ev.strand})])
reg = reg[reg.e > reg.s]
reg.to_csv(f"{U}/_regions.bed", sep="\t", header=False, index=False)

man = pd.read_csv(f"{PB}/01_encode_metadata/download_manifest.tsv", sep="\t", dtype=str)
pk = man[(man.assay == "eCLIP") & (man.bio_reps == "1;2") & (man.cell_line.isin(["K562", "HepG2"]))]
print("IDR-merged peak files:", len(pk), "RBPs:", pk.target.nunique(), file=sys.stderr)
with open(f"{U}/_peaks.bed", "w") as out:
    for r in pk.itertuples():
        for line in gzip.open(f"{PB}/02_encode_download/{r.dest}", "rt"):
            f = line.rstrip("\n").split("\t")
            out.write(f"{f[0]}\t{f[1]}\t{f[2]}\t{r.target}|{r.cell_line}\t0\t{f[5]}\t{f[6]}\t{f[7]}\n")
cmd = (f"sort -k1,1 -k2,2n {U}/_regions.bed > {U}/_regions.sorted.bed && sort -k1,1 -k2,2n {U}/_peaks.bed > {U}/_peaks.sorted.bed && "
       f"bedtools intersect -a {U}/_regions.sorted.bed -b {U}/_peaks.sorted.bed -s -wa -wb -sorted | cut -f4,10,13,14 > {U}/_hits.tsv && "
       f"rm {U}/_regions.bed {U}/_regions.sorted.bed {U}/_peaks.bed {U}/_peaks.sorted.bed")
subprocess.run(["bash", "-c", "set -o pipefail; " + cmd], check=True)
h = pd.read_csv(f"{U}/_hits.tsv", sep="\t", header=None, names=["name", "exp", "sig", "nlp"])
h[["event_id", "region"]] = h.name.str.split("|", expand=True); h[["rbp", "cell_line"]] = h.exp.str.split("|", expand=True)
F = h.pivot_table(index=["event_id", "rbp", "cell_line"], columns="region", values=["sig", "nlp"], aggfunc="max", fill_value=0)
F.columns = [f"{a}_{b}" for a, b in F.columns]; F = F.reset_index()
F.to_parquet(f"{U}/eclip_event_features.parquet", index=False)
subprocess.run(["rm", f"{U}/_hits.tsv"])
print("rows:", len(F), "| events with any peak:", F.event_id.nunique(), "/", len(ev), "| RBP x cell line:", F.groupby(["rbp", "cell_line"]).ngroups, file=sys.stderr)
print(F.groupby(["rbp", "cell_line"]).size().describe().round(0).to_string(), file=sys.stderr)
print((F[["sig_up", "sig_seg", "sig_dn"]] > 0).mean().round(3).to_string(), file=sys.stderr)
