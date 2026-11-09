"""Sequence / conservation features for every event in events.parquet (coordinates are GRCh38, 0-based half-open as rMATS gives them).
The "alternative segment" is: SE = the cassette exon; A3SS/A5SS = long exon minus short exon; MXE = 1st exon; RI = the retained intron.
Features: segment length, GC, upstream/downstream intron lengths, MaxEntScan 5' (9-mer) and 3' (23-mer) scores at the segment's
donor and acceptor (strand-aware), phastCons100way mean over the segment and over 100 nt of each flanking intron.
A3SS/A5SS: only the alternative splice site of the long exon is scored (A3SS acceptor, A5SS donor); the other score and the
intron on that side are NaN; the single intron length is the gap to the flanking exon. RI: donor/acceptor are the intron's own
5'/3' splice sites; flanking introns and flank conservation are NaN.
Usage: python 02_event_features.py <universe_dir>   -> <universe_dir>/event_features.parquet"""
import sys, os, subprocess, numpy as np, pandas as pd, pyBigWig
from pyfaidx import Fasta
PB = "/well/blagden/users/eic789/PerturbBridge"; U = sys.argv[1]
FA = Fasta("/well/blagden/users/eic789/LncRNASPLICE/GRCh38.v47.primary_assembly.genome.fa")
BW = pyBigWig.open(f"{PB}/06_features/resources/hg38.phastCons100way.bw")
MAXENT = f"{PB}/06_features/resources/maxentpy_src/maxentpy-master/maxent"
ev = pd.read_parquet(f"{U}/events.parquet")
def segs(r):
    c = [int(v) for v in r.coords.split(",")]
    if r.type == "SE": s, e, upE, dnS = c[0], c[1], c[3], c[4]
    elif r.type == "MXE": s, e, upE, dnS = c[0], c[1], c[5], c[6]
    elif r.type == "RI": s, e, upE, dnS = c[3], c[4], c[3], c[4]          # intron between upstreamEE and downstreamES
    else:  # A3SS/A5SS: long minus short exon; flanking exon on the other side
        L0, L1, S0, S1, F0, F1 = c
        if L0 == S0: s, e = S1, L1          # alternative part at the 3' end (genomic) of the long exon
        else: s, e = L0, S0
        upE, dnS = (F1, s) if F1 <= s else (e, F0)
    return s, e, upE, dnS
def rc(x): return x.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
def seq(ch, s, e): return FA[ch][max(0, s):e].seq.upper()
def bw_mean(ch, s, e):
    if ch not in BW.chroms() or e <= s: return np.nan
    v = BW.stats(ch, max(0, s), min(e, BW.chroms()[ch]), type="mean")[0]; return np.nan if v is None else v
rows, don, acc = [], [], []
for i, r in enumerate(ev.itertuples()):
    s, e, upE, dnS = segs(r); ch, st = r.chr, r.strand
    if ch not in FA.keys(): rows.append((r.event_id, s, e, e - s, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)); don.append("N"*9); acc.append("N"*23); continue
    sg = seq(ch, s, e); gc = (sg.count("G") + sg.count("C")) / max(1, len(sg))
    if r.type == "RI":  # segment is the intron: donor at its 5' end, acceptor at its 3' end (transcript orientation)
        if st == "+": d = seq(ch, s - 3, s + 6); a = seq(ch, e - 20, e + 3)
        else: d = rc(seq(ch, e - 6, e + 3)); a = rc(seq(ch, s - 3, s + 20))
        upI = dnI = np.nan; pc_up, pc_dn = np.nan, np.nan
    elif st == "+":
        d = seq(ch, e - 3, e + 6); a = seq(ch, s - 20, s + 3); upI, dnI = s - upE, dnS - e
        pc_up, pc_dn = bw_mean(ch, s - 100, s), bw_mean(ch, e, e + 100)
    else:
        d = rc(seq(ch, s - 6, s + 3)); a = rc(seq(ch, e - 3, e + 20)); upI, dnI = dnS - e, s - upE
        pc_up, pc_dn = bw_mean(ch, e, e + 100), bw_mean(ch, s - 100, s)
    if r.type in ("A3SS", "A5SS"):  # only one side of the alternative segment is a splice site; the intron is the gap to the flanking exon
        L0, L1, S0, S1, F0, F1 = [int(v) for v in r.coords.split(",")]; gap = L0 - F1 if F1 <= L0 else F0 - L1
        if r.type == "A3SS": upI, dnI, d = gap, np.nan, "N" * 9
        else: upI, dnI, a = np.nan, gap, "N" * 23
    don.append(d if len(d) == 9 else "N"*9); acc.append(a if len(a) == 23 else "N"*23)
    rows.append((r.event_id, s, e, e - s, gc, upI, dnI, bw_mean(ch, s, e), pc_up, pc_dn))
    if i % 100000 == 0: print(i, file=sys.stderr, flush=True)
F = pd.DataFrame(rows, columns=["event_id", "seg_start", "seg_end", "seg_len", "seg_gc", "up_intron_len", "down_intron_len", "phastcons_seg", "phastcons_up100", "phastcons_down100"])
def maxent(seqs, script):  # MaxEntScan cannot score sequences containing N: score the valid ones, NaN for the rest
    ok = [i for i, q in enumerate(seqs) if "N" not in q]; inp = f"{U}/_me_in.txt"; open(inp, "w").write("\n".join(seqs[i] for i in ok) + "\n")
    o = subprocess.run(f"cd {MAXENT} && perl {script} {os.path.abspath(inp)}", shell=True, capture_output=True, text=True, check=True).stdout
    sc = [float(l.split("\t")[1]) for l in o.strip().split("\n")]; os.remove(inp); assert len(sc) == len(ok)
    out = np.full(len(seqs), np.nan); out[ok] = sc; return out
F["maxent_donor"] = maxent(don, "score5.pl"); F["maxent_acceptor"] = maxent(acc, "score3.pl")
F["donor_9mer"] = don; F["acceptor_23mer"] = acc
F.to_parquet(f"{U}/event_features.parquet", index=False)
print(F.drop(columns=["event_id", "donor_9mer", "acceptor_23mer"]).describe().round(2).T.to_string())
