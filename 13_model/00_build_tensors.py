"""Pack everything the forward model needs for one cell line into a single npz (design 3):
  events   : the cell line's training events (train_table_<cl>.parquet), event features (type one-hot, log lengths, GC, MaxEntScan,
             phastCons, missing flags) and the length-normalised 5-mer spectrum of up/seg/dn (3 x 1024)
  rbps     : ESM2 (1280), GenePT (1536), motif 5-mer vector (1024, zero + flag if absent) per target of this cell line
  eclip    : per (target with eCLIP in this cell line) x event: sig_up/seg/dn, nlp_up/seg/dn (zeros where no peak)
  rows     : one per (event, run): event index, run index, psi_ctrl, y = control-centred dPSI (dpsi_cc), raw dpsi, sig flag
  runs     : run -> target index, assay (CRISPR=1), tier, fold, control
Usage: python 00_build_tensors.py <cell_line>  -> 13_model/tensors_<cell_line>.npz"""
import sys, pickle, numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; U = f"{PB}/10_event_universe/encode_full"; cl = sys.argv[1]
t = pd.read_parquet(f"{U}/train_table_{cl}.parquet", columns=["event_id", "run", "target", "assay", "psi_ctrl", "dpsi", "dpsi_cc", "fdr"])
t = t[np.isfinite(t.dpsi_cc)]
events = np.array(sorted(t.event_id.unique())); epos = {e: i for i, e in enumerate(events)}
tiers = pd.read_csv(f"{PB}/11_splits/run_tiers.tsv", sep="\t").set_index("run"); runs = np.array(sorted(t.run.unique())); rpos = {r: i for i, r in enumerate(runs)}
targets = np.array(sorted(t.target.unique())); tpos = {x: i for i, x in enumerate(targets)}
# event features
ev = pd.read_parquet(f"{U}/events.parquet", columns=["event_id", "type"]).merge(pd.read_parquet(f"{U}/event_features.parquet"), on="event_id").set_index("event_id").loc[events]
types = ["SE", "A3SS", "A5SS", "MXE", "RI"]; X = [np.stack([(ev.type == x).values for x in types], 1).astype(np.float32)]
def col(v, log=False):
    v = v.values.astype(np.float64); m = np.isfinite(v); w = np.where(m, v, 0.0)
    if log: w = np.log10(w + 1)
    return np.stack([w, (~m).astype(float)], 1).astype(np.float32)
for c, lg in [("seg_len", True), ("seg_gc", False), ("up_intron_len", True), ("down_intron_len", True), ("maxent_donor", False), ("maxent_acceptor", False), ("phastcons_seg", False), ("phastcons_up100", False), ("phastcons_down100", False)]:
    X.append(col(ev[c], lg))
Xe = np.hstack(X); names = types + [f"{c}{s}" for c in ["seg_len", "seg_gc", "up_intron_len", "down_intron_len", "maxent_donor", "maxent_acceptor", "phastcons_seg", "phastcons_up100", "phastcons_down100"] for s in ("", "_missing")]
z = np.load(f"{U}/kmer5_counts.npz"); kpos = {e: i for i, e in enumerate(z["event_id"])}; kidx = np.array([kpos[e] for e in events])
K = np.stack([z[r][kidx].astype(np.float32) for r in ("up", "seg", "dn")], 1); K = K / (K.sum(2, keepdims=True) + 1e-6)   # events x 3 x 1024, length-normalised
# RBP features
ix = pd.read_csv(f"{PB}/06_features/rbp_annotation/rbp_esm2_650M_index.tsv", sep="\t"); E = np.load(f"{PB}/06_features/rbp_annotation/rbp_esm2_650M_mean.npy"); esm = dict(zip(ix.target, E))
gp = pickle.load(open(f"{PB}/06_features/rbp_annotation/genept/GenePT_emebdding_v2/GenePT_gene_embedding_ada_text.pickle", "rb"))
mv = pd.read_csv(f"{PB}/06_features/rbp_annotation/rbp_motif_5mer.tsv", sep="\t", index_col=0).drop(columns="source")
R_esm = np.array([esm[x] for x in targets], dtype=np.float32); R_gp = np.array([np.asarray(gp[x]) if x in gp else np.zeros(1536) for x in targets], dtype=np.float32)
R_mot = np.array([mv.loc[x].values if x in mv.index else np.zeros(1024) for x in targets], dtype=np.float32); R_has_mot = np.array([x in mv.index for x in targets], dtype=np.float32)
# eCLIP per (target, event) for targets with eCLIP in this cell line
ec = pd.read_parquet(f"{U}/eclip_event_features.parquet"); ec = ec[(ec.cell_line == cl) & ec.event_id.isin(epos) & ec.rbp.isin(tpos)]
ec_targets = np.array(sorted(ec.rbp.unique())); ecpos = {x: i for i, x in enumerate(ec_targets)}
EC = np.zeros((len(ec_targets), len(events), 6), dtype=np.float32)
EC[[ecpos[x] for x in ec.rbp], [epos[e] for e in ec.event_id]] = ec[["sig_up", "sig_seg", "sig_dn", "nlp_up", "nlp_seg", "nlp_dn"]].values
EC[..., 3:] = np.log1p(EC[..., 3:])
target_ecidx = np.array([ecpos.get(x, -1) for x in targets])
# rows and runs
rows = pd.DataFrame({"e": t.event_id.map(epos).values, "r": t.run.map(rpos).values, "psi_ctrl": t.psi_ctrl.values, "y": t.dpsi_cc.values, "y_raw": t.dpsi.values,
                     "sig": ((t.fdr < 0.05) & (t.dpsi.abs() > 0.1)).values})
ri = tiers.loc[runs]; run_target = np.array([tpos[x] for x in ri.target]); run_assay = ri.assay.str.startswith("CRISPR").values.astype(np.float32)
run_tier = (ri.tier == "A").values.astype(np.float32); run_fold = ri.fold.values; run_control = ri.control.values
np.savez_compressed(f"{PB}/13_model/tensors_{cl}.npz", events=events, Xe=Xe, Xe_names=np.array(names), K=K, targets=targets, R_esm=R_esm, R_gp=R_gp, R_mot=R_mot, R_has_mot=R_has_mot,
                    ec_targets=ec_targets, EC=EC, target_ecidx=target_ecidx, runs=runs, run_target=run_target, run_assay=run_assay, run_tier=run_tier, run_fold=run_fold, run_control=run_control,
                    row_e=rows.e.values.astype(np.int32), row_r=rows.r.values.astype(np.int32), row_psi=rows.psi_ctrl.values.astype(np.float32), row_y=rows.y.values.astype(np.float32),
                    row_yraw=rows.y_raw.values.astype(np.float32), row_sig=rows.sig.values)
print(cl, "events", len(events), "runs", len(runs), "targets", len(targets), "rows", len(rows), "| Xe", Xe.shape, "K", K.shape, "EC", EC.shape, "| targets with eCLIP", len(ec_targets), "with motif", int(R_has_mot.sum()))
