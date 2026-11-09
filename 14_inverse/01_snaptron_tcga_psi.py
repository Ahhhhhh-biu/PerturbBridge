"""Per-sample PSI for the model's event universe in TCGA, from Snaptron (tcgav2 compilation, recount2 junction counts).
For each event (SE, A3SS, A5SS, MXE; RI skipped because it needs intron-body coverage) the defining junctions are queried with
exact=1 (rMATS 0-based exon starts -> 1-based intron coordinates: intron = [exonEnd+1, nextExonStart_0base]); counts are kept
only for samples of the projects in PROJECTS. PSI per sample = mean inclusion-junction count / (mean inclusion + skipping),
reported with the total junction coverage so callers can filter (>= 10 reads).
Usage: python 01_snaptron_tcga_psi.py <task_id> <n_tasks>  -> 14_inverse/tcga_psi/part_<task>.parquet (event_id, rail_id, inc, skip, psi)
Bulk POST queries (10 events = 30-40 regions per request, contains=1, exact coordinates filtered locally); 2 concurrent streams per shard; resumable (part_<task>.done_events.txt). Retries up to 10 times."""
import sys, time, urllib.request, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
PB = "/well/blagden/users/eic789/PerturbBridge"; O = f"{PB}/14_inverse/tcga_psi"; task, ntask = int(sys.argv[1]), int(sys.argv[2])
PROJECTS = {"TCGA-UVM", "TCGA-LUAD", "TCGA-BRCA", "TCGA-OV", "TCGA-LAML", "TCGA-SKCM"}
S = pd.read_csv(f"{PB}/14_inverse/tcga_samples.tsv", sep="\t", dtype=str, usecols=["rail_id", "gdc_cases.project.project_id"]); keep = set(S.loc[S["gdc_cases.project.project_id"].isin(PROJECTS), "rail_id"].astype(int))
ev = pd.read_parquet(f"{PB}/10_event_universe/encode_full/events.parquet")
E = sorted(set(np.load(f"{PB}/13_model/tensors_HepG2.npz", allow_pickle=True)["events"]) | set(np.load(f"{PB}/13_model/tensors_K562.npz", allow_pickle=True)["events"]))
ev = ev[ev.event_id.isin(E) & (ev.type != "RI")].sort_values("event_id").reset_index(drop=True); ev = ev.iloc[task::ntask]
def junctions(r):
    c = list(map(int, r.coords.split(","))); p = r.strand == "+"
    if r.type == "SE": eS, eE, uS, uE, dS, dE = c; return [(uE + 1, eS), (eE + 1, dS)], [(uE + 1, dS)]
    if r.type == "A3SS": lS, lE, sS, sE, fS, fE = c; return ([(fE + 1, lS)], [(fE + 1, sS)]) if p else ([(lE + 1, fS)], [(sE + 1, fS)])
    if r.type == "A5SS": lS, lE, sS, sE, fS, fE = c; return ([(lE + 1, fS)], [(sE + 1, fS)]) if p else ([(fE + 1, lS)], [(fE + 1, sS)])
    if r.type == "MXE": aS, aE, bS, bE, uS, uE, dS, dE = c; return [(uE + 1, aS), (aE + 1, dS)], [(uE + 1, bS), (bE + 1, dS)]
import json, urllib.parse
def fetch_bulk(regions):
    """One POST with all regions (contains=1); returns {(chr,start,end): {rail_id: count}} for junctions whose coordinates match a requested region exactly."""
    want = set(regions); fields = [dict(intervals=[f"{c}:{s}-{e}"], contains=1) for c, s, e in regions]
    body = urllib.parse.urlencode({"fields": json.dumps(fields)}).encode()
    for k in range(10):
        try:
            txt = urllib.request.urlopen(urllib.request.Request("https://snaptron.cs.jhu.edu/tcgav2/snaptron", data=body), timeout=600).read().decode()
            out = {}
            for line in txt.split("\n"):
                if not line.startswith("tcgav2"): continue
                f = line.split("\t"); key = (f[2], int(f[3]), int(f[4]))
                if key not in want or key in out: continue
                d = {}
                for x in f[12].split(",")[1:]:
                    rid, cnt = x.split(":"); rid = int(rid)
                    if rid in keep: d[rid] = d.get(rid, 0) + int(cnt)
                out[key] = d
            return out
        except Exception as ex:
            time.sleep(min(90, 3 * 2 ** k)); err = ex
    raise RuntimeError(f"bulk query failed ({len(regions)} regions): {err}")
def chunk_rows(evs):
    """evs: list of event rows -> PSI rows for all of them from one bulk query."""
    J = {r.event_id: junctions(r) for r in evs}; regions = sorted({(r.chr, s, e) for r in evs for grp in J[r.event_id] for s, e in grp})
    got = fetch_bulk(regions); rows = []
    for r in evs:
        inc, skp = J[r.event_id]; I = [got.get((r.chr, s, e), {}) for s, e in inc]; K = [got.get((r.chr, s, e), {}) for s, e in skp]
        for rid in set().union(*I, *K):
            a = np.mean([d.get(rid, 0) for d in I]); b = np.mean([d.get(rid, 0) for d in K])
            if a + b >= 5: rows.append((r.event_id, rid, a, b, a / (a + b)))
    return rows
import os
t0 = time.time(); out = f"{O}/part_{task:03d}.parquet"; done_file = f"{O}/part_{task:03d}.done_events.txt"
done = set(open(done_file).read().split()) if os.path.exists(done_file) else set(); parts = [pd.read_parquet(out)] if os.path.exists(out) and done else []
todo = ev[~ev.event_id.isin(done)]; print(f"task {task}: {len(ev)} events, {len(done)} already done, {len(todo)} to do", flush=True)
def flush(res, finished):
    global parts
    df = pd.DataFrame(res, columns=["event_id", "rail_id", "inc", "skip", "psi"]); df["rail_id"] = df.rail_id.astype(np.int32); df["inc"] = df.inc.astype(np.float32); df["skip"] = df.skip.astype(np.float32); df["psi"] = df.psi.astype(np.float32)
    parts.append(df); pd.concat(parts, ignore_index=True).to_parquet(out, index=False); parts = [pd.read_parquet(out)]
    with open(done_file, "a") as f: f.write("\n".join(finished) + "\n")
res, finished = [], []; CH = 10; todo_rows = list(todo.itertuples()); chunks = [todo_rows[i:i + CH] for i in range(0, len(todo_rows), CH)]
with ThreadPoolExecutor(2) as ex:
    for i, (ch, rows) in enumerate(zip(chunks, ex.map(chunk_rows, chunks))):
        res.extend(rows); finished += [r.event_id for r in ch]
        if (i + 1) % 10 == 0: flush(res, finished); res, finished = [], []; print(f"task {task}: {(i + 1) * CH}/{len(todo)} events, {time.time() - t0:.0f}s", flush=True)
if finished: flush(res, finished)
print("DONE", task, len(ev), "events", f"{time.time() - t0:.0f}s")
