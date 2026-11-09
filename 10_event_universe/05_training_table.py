"""Training delta-PSI table per cell line: events covered (both groups >=10 reads) in >= COVFRAC of that cell line's runs,
one row per (event, run) that is covered in that run. dpsi = PSI(KD) - PSI(control) as rMATS gives it (b1 = KD).
Usage: python 05_training_table.py <universe_dir>  -> <universe_dir>/train_table_<cell_line>.parquet"""
import sys, pandas as pd, pyarrow.parquet as pq
PB = "/well/blagden/users/eic789/PerturbBridge"; U = sys.argv[1]; MINCOV = 10; COVFRAC = 0.5
pairs = pd.read_csv(f"{PB}/09_bam_layer/kd_control_pairs_all.tsv", sep="\t"); pairs["run"] = pairs.kd + "_vs_" + pairs.control
cov = pd.read_csv(f"{U}/summary_event_coverage.tsv", sep="\t")
t = pq.read_table(f"{U}/per_run.parquet", columns=["event_id", "run", "cov1", "cov2", "psi1", "psi2", "dpsi", "fdr", "psi1_reps", "psi2_reps"], read_dictionary=["event_id", "run"]).to_pandas()
t = t[(t.cov1 >= MINCOV) & (t.cov2 >= MINCOV)]; t["run"] = t.run.astype(str); t["event_id"] = t.event_id.astype(str)
t = t.merge(pairs[["run", "cell_line", "target", "assay"]], on="run")
for cl in ["HepG2", "K562"]:
    keep = set(cov.loc[cov[f"covfrac_{cl}"] >= COVFRAC, "event_id"]); d = t[(t.cell_line == cl) & t.event_id.isin(keep)]
    d = d.rename(columns={"psi1": "psi_kd", "psi2": "psi_ctrl", "psi1_reps": "psi_kd_reps", "psi2_reps": "psi_ctrl_reps"})
    d.to_parquet(f"{U}/train_table_{cl}.parquet", index=False)
    print(cl, "events:", len(keep), "runs:", d.run.nunique(), "targets:", d.target.nunique(), "rows:", len(d),
          "| rows per event median:", int(d.groupby("event_id").size().median()), "| |dpsi|>0.1:", round((d.dpsi.abs() > 0.1).mean(), 4), "| sig:", round(((d.fdr < 0.05) & (d.dpsi.abs() > 0.1)).mean(), 4))
