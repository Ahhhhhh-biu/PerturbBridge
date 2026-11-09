"""Build a coordinate-keyed event universe from many rMATS-turbo runs and a long per-run table.
rMATS event IDs are assigned per run, so events are keyed here by (type, chr, strand, exon coordinates) instead.
Usage: python 01_build_event_universe.py <out_dir> <rmats_dir> [<rmats_dir> ...]   (each dir holds *.MATS.JC.txt, b1.txt, b2.txt)
Outputs: <out_dir>/events.parquet (one row per distinct event), <out_dir>/per_run.parquet (one row per event x run)."""
import sys, os, numpy as np, pandas as pd
out, dirs = sys.argv[1], sys.argv[2:]
COORD = {"SE": ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
         "A3SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
         "A5SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
         "MXE": ["1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
         "RI": ["riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"]}
def rep_sum(col): return col.str.split(",").apply(lambda v: sum(int(x) for x in v))
def rep_mean_na(col): return col.str.split(",").apply(lambda v: np.nanmean([float(x) if x != "NA" else np.nan for x in v]) if any(x != "NA" for x in v) else np.nan)
events, runs = {}, []
for d in dirs:
    run = os.path.basename(d.rstrip("/"))
    for t, cols in COORD.items():
        f = f"{d}/{t}.MATS.JC.txt"
        if not os.path.exists(f): continue
        x = pd.read_csv(f, sep="\t", dtype={"IJC_SAMPLE_1": str, "SJC_SAMPLE_1": str, "IJC_SAMPLE_2": str, "SJC_SAMPLE_2": str, "IncLevel1": str, "IncLevel2": str})
        key = t + "|" + x["chr"] + "|" + x["strand"] + "|" + x[cols].astype(str).agg("|".join, axis=1)
        for k, g, s in zip(key, x.GeneID, x.geneSymbol):
            if k not in events: events[k] = (g, s)
        n1 = x.IJC_SAMPLE_1.str.count(",") + 1; n2 = x.IJC_SAMPLE_2.str.count(",") + 1
        runs.append(pd.DataFrame({"key": key, "run": run,
            "ijc1": rep_sum(x.IJC_SAMPLE_1), "sjc1": rep_sum(x.SJC_SAMPLE_1), "ijc2": rep_sum(x.IJC_SAMPLE_2), "sjc2": rep_sum(x.SJC_SAMPLE_2),
            "cov1": (rep_sum(x.IJC_SAMPLE_1) + rep_sum(x.SJC_SAMPLE_1)) / n1, "cov2": (rep_sum(x.IJC_SAMPLE_2) + rep_sum(x.SJC_SAMPLE_2)) / n2,
            "psi1": rep_mean_na(x.IncLevel1), "psi2": rep_mean_na(x.IncLevel2), "dpsi": x.IncLevelDifference.astype(float),
            "pvalue": x.PValue.astype(float), "fdr": x.FDR.astype(float), "psi1_reps": x.IncLevel1, "psi2_reps": x.IncLevel2}))
    print("loaded", run, file=sys.stderr)
ev = pd.DataFrame([(k, *k.split("|")[:3], k.split("|")[3:], g, s) for k, (g, s) in events.items()], columns=["key", "type", "chr", "strand", "coords", "gene_id", "gene_symbol"])
ev["coords"] = ev.coords.apply(lambda v: ",".join(v)); ev["_pos"] = ev.coords.str.split(",").str[0].astype(int)
ev = ev.sort_values(["type", "chr", "_pos", "coords"]).drop(columns="_pos").reset_index(drop=True)
ev["event_id"] = ev.type + "_" + ev.groupby("type").cumcount().astype(str).str.zfill(6)
pr = pd.concat(runs, ignore_index=True).merge(ev[["key", "event_id"]], on="key").drop(columns="key")
os.makedirs(out, exist_ok=True); ev.to_parquet(f"{out}/events.parquet", index=False); pr.to_parquet(f"{out}/per_run.parquet", index=False)
print("events:", len(ev), ev.type.value_counts().to_dict()); print("runs:", pr.run.nunique(), "rows:", len(pr))
print("events seen in all runs:", (pr.groupby("event_id").run.nunique() == pr.run.nunique()).sum(), "| seen in >=half:", (pr.groupby("event_id").run.nunique() >= pr.run.nunique() / 2).sum())
