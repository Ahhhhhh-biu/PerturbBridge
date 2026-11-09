#!/usr/bin/env python
"""Event-level eCLIP features: for every SUPPA event, windows of ±300 nt around each splice-site coordinate
encoded in the event id; intersect with each RBP's IDR-merged eCLIP narrowPeak (per cell line).
Output: 06_features/eclip/eclip_event_features.parquet  (event_id, rbp, cell_line, max_signal, max_neglog10p, n_sites_hit)"""
import pandas as pd, numpy as np, re, subprocess, os, gzip, sys
PB = "/well/blagden/users/eic789/PerturbBridge"; OUT = f"{PB}/06_features/eclip"; W = 300
ioe = pd.read_csv(f"{PB}/03_suppa_psi/events/gencode_v29_all_strict.ioe", sep="\t")
rows = []
for eid in ioe.event_id:
    body = eid.split(";", 1)[1]                      # e.g. SE:chr1:100-200:300-400:+
    parts = body.split(":"); chrom = parts[1]; strand = parts[-1]
    coords = sorted({int(x) for seg in parts[2:-1] for x in re.findall(r"\d+", seg)})
    for k, c in enumerate(coords):
        rows.append((chrom, max(0, c - W), c + W, eid, k, strand))
win = pd.DataFrame(rows, columns=["chrom", "start", "end", "event_id", "site", "strand"])
win.to_csv(f"{OUT}/event_site_windows.bed", sep="\t", header=False, index=False)
subprocess.run(f"sort -k1,1 -k2,2n {OUT}/event_site_windows.bed > {OUT}/event_site_windows.sorted.bed", shell=True, check=True)
print("windows:", len(win), "events:", win.event_id.nunique(), file=sys.stderr)

man = pd.read_csv(f"{PB}/01_encode_metadata/download_manifest.tsv", sep="\t", dtype=str)
pk = man[(man.assay == "eCLIP") & (man.bio_reps == "1;2") & (man.cell_line.isin(["K562", "HepG2"]))]
print("IDR-merged peak files:", len(pk), "RBPs:", pk.target.nunique(), file=sys.stderr)
feats = []
for r in pk.itertuples():
    bed = f"{PB}/02_encode_download/{r.dest}"
    cmd = f"bedtools intersect -a {OUT}/event_site_windows.sorted.bed -b <(zcat {bed} | cut -f1-3,7,8) -wa -wb -sorted 2>/dev/null || bedtools intersect -a {OUT}/event_site_windows.sorted.bed -b <(zcat {bed} | cut -f1-3,7,8) -wa -wb"
    res = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if not res.stdout.strip(): continue
    d = pd.read_csv(pd.io.common.StringIO(res.stdout), sep="\t", header=None,
                    usecols=[3, 4, 9, 10], names=["event_id", "site", "signal", "neglog10p"])
    g = d.groupby("event_id").agg(max_signal=("signal", "max"), max_neglog10p=("neglog10p", "max"), n_sites_hit=("site", "nunique")).reset_index()
    g["rbp"] = r.target; g["cell_line"] = r.cell_line; g["file"] = r.file
    feats.append(g)
F = pd.concat(feats, ignore_index=True)
F.to_parquet(f"{OUT}/eclip_event_features.parquet", index=False)
print("feature rows:", len(F), "| RBP×cell_line combos:", F.groupby(["rbp", "cell_line"]).ngroups, file=sys.stderr)
print(F.groupby(["rbp", "cell_line"]).size().describe().to_string(), file=sys.stderr)
