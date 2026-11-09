"""Raw-sequence windows for the 1D-CNN event encoder: 300 nt upstream intron | alternative segment (first 150 + last 150 nt; shorter
segments centred and N-padded to 300) | 300 nt downstream intron, in transcript orientation (reverse-complemented on the minus strand).
Encoded as uint8 codes A0 C1 G2 T3 N4, shape events x 900, row order = 13_model/tensors_<cell>.npz['events'].
Usage: python 00_event_sequences.py <cell_line> -> 13_model/tensors_seq_<cell_line>.npz"""
import sys, numpy as np, pandas as pd
from pyfaidx import Fasta
PB = "/well/blagden/users/eic789/PerturbBridge"; U = f"{PB}/10_event_universe/encode_full"; cl = sys.argv[1]; W = 300; H = 150
FA = Fasta("/well/blagden/users/eic789/LncRNASPLICE/GRCh38.v47.primary_assembly.genome.fa")
events = np.load(f"{PB}/13_model/tensors_{cl}.npz", allow_pickle=True)["events"]
ev = pd.read_parquet(f"{U}/events.parquet", columns=["event_id", "chr", "strand"]).merge(pd.read_parquet(f"{U}/event_features.parquet", columns=["event_id", "seg_start", "seg_end"]), on="event_id").set_index("event_id").loc[events]
LUT = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate("ACGT"): LUT[ord(b)] = i; LUT[ord(b.lower())] = i
COMP = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
def fetch(chrom, s, e):   # 0-based half-open, clipped, N-padded
    L = len(FA[chrom]); a, b = max(0, s), min(L, e); seq = LUT[np.frombuffer(str(FA[chrom][a:b]).encode(), dtype=np.uint8)] if b > a else np.zeros(0, np.uint8)
    return np.concatenate([np.full(a - s, 4, np.uint8), seq, np.full(e - b, 4, np.uint8)])
out = np.zeros((len(ev), 3 * W), dtype=np.uint8)
for i, r in enumerate(ev.itertuples()):
    s0, s1 = int(r.seg_start), int(r.seg_end); seglen = s1 - s0
    if seglen >= 2 * H: seg = np.concatenate([fetch(r.chr, s0, s0 + H), fetch(r.chr, s1 - H, s1)])
    else: pad = 2 * H - seglen; seg = np.concatenate([np.full(pad // 2, 4, np.uint8), fetch(r.chr, s0, s1), np.full(pad - pad // 2, 4, np.uint8)])
    left = fetch(r.chr, s0 - W, s0); right = fetch(r.chr, s1, s1 + W); x = np.concatenate([left, seg, right])
    if r.strand == "-": x = COMP[x[::-1]]
    out[i] = x
    if i % 10000 == 0: print(i, flush=True)
np.savez_compressed(f"{PB}/13_model/tensors_seq_{cl}.npz", events=events, seq=out); print("DONE", cl, out.shape, "N fraction", float((out == 4).mean()))
