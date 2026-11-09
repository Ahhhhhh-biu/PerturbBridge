"""5-mer counts for the three event regions used by the eCLIP features (seg = alternative segment, up = 300 nt upstream,
dn = 300 nt downstream, transcript orientation; reverse-complemented on the minus strand). Counts are over the RNA sense
strand so they can be dotted with RBNS 5-mer enrichment vectors (06_features/rbns/rbns_5mer_R.tsv) to give per-event,
per-region motif scores for any RBP with RBNS data, and used directly as the k-mer spectrum of the event encoder.
Output: <universe_dir>/kmer5_counts.npz with arrays up, seg, dn (uint16, events x 1024, row order = events.parquet), kmers, event_id.
Usage: python 06_kmer_counts.py <universe_dir>"""
import sys, itertools, numpy as np, pandas as pd
from pyfaidx import Fasta
U = sys.argv[1]; W = 300; K = 5
FA = Fasta("/well/blagden/users/eic789/LncRNASPLICE/GRCh38.v47.primary_assembly.genome.fa")
ev = pd.read_parquet(f"{U}/events.parquet", columns=["event_id", "chr", "strand"]).merge(pd.read_parquet(f"{U}/event_features.parquet", columns=["event_id", "seg_start", "seg_end"]), on="event_id")
LUT = np.full(256, 4, dtype=np.int64); 
for i, b in enumerate("ACGT"): LUT[ord(b)] = i; LUT[ord(b.lower())] = i
kmers = ["".join(p) for p in itertools.product("ACGT", repeat=K)]
def count(seq, rc):
    a = LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]
    if rc: a = (3 - a)[::-1]; a[a < 0] = 4
    if len(a) < K: return np.zeros(4 ** K, dtype=np.uint16)
    win = np.lib.stride_tricks.sliding_window_view(a, K); ok = (win < 4).all(1)
    code = (win[ok] * (4 ** np.arange(K - 1, -1, -1))).sum(1); return np.bincount(code, minlength=4 ** K).astype(np.uint16)
out = {r: np.zeros((len(ev), 4 ** K), dtype=np.uint16) for r in ("up", "seg", "dn")}
for i, r in enumerate(ev.itertuples()):
    if r.chr not in FA.keys(): continue
    L = len(FA[r.chr]); s, e = int(r.seg_start), int(r.seg_end); rc = r.strand == "-"
    seg = FA[r.chr][max(0, s):e].seq; left = FA[r.chr][max(0, s - W):s].seq; right = FA[r.chr][e:min(L, e + W)].seq
    out["seg"][i] = count(seg, rc); out["up" if not rc else "dn"][i] = count(left, rc); out["dn" if not rc else "up"][i] = count(right, rc)
    if i % 100000 == 0: print(i, file=sys.stderr, flush=True)
np.savez_compressed(f"{U}/kmer5_counts.npz", kmers=np.array(kmers, dtype=str), event_id=np.array(list(ev.event_id), dtype=str), **out)
print("done", {k: (v.shape, int(v.sum(1).mean())) for k, v in out.items()})
