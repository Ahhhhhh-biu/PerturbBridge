"""Coverage and effect-size summary of the full ENCODE event universe.
Usage: python 04_universe_summary.py <universe_dir>  -> <universe_dir>/summary_per_run.tsv, summary_event_coverage.tsv (+ stdout)"""
import sys, numpy as np, pandas as pd, pyarrow.parquet as pq
PB = "/well/blagden/users/eic789/PerturbBridge"; U = sys.argv[1]; MINCOV = 10
pairs = pd.read_csv(f"{PB}/09_bam_layer/kd_control_pairs_all.tsv", sep="\t"); pairs["run"] = pairs.kd + "_vs_" + pairs.control
t = pq.read_table(f"{U}/per_run.parquet", columns=["event_id", "run", "cov1", "cov2", "dpsi", "fdr"], read_dictionary=["event_id", "run"]).to_pandas()
ev = pd.read_parquet(f"{U}/events.parquet", columns=["event_id", "type"])
t["covered"] = (t.cov1 >= MINCOV) & (t.cov2 >= MINCOV); t["sig"] = t["covered"] & (t.fdr < 0.05) & (t.dpsi.abs() > 0.1)
t["run"] = t.run.astype(str); t = t.merge(pairs[["run", "cell_line", "target", "assay"]], on="run")
pr = t.groupby(["run", "cell_line", "target", "assay"], observed=True).agg(n_events=("event_id", "size"), n_cov=("covered", "sum"), n_sig=("sig", "sum")).reset_index()
pr["sig_frac"] = (pr.n_sig / pr.n_cov).round(4); pr.sort_values("n_sig", ascending=False).to_csv(f"{U}/summary_per_run.tsv", sep="\t", index=False)
print("== runs per cell line / assay ==\n", pr.groupby(["cell_line", "assay"]).size().to_string())
print("\n== per-run covered events and significant events (FDR<0.05, |dPSI|>0.1) ==")
print(pr.groupby(["cell_line", "assay"])[["n_cov", "n_sig", "sig_frac"]].describe(percentiles=[.1, .5, .9]).round(3).T.to_string())
print("\n== top 15 runs by n_sig ==\n", pr.sort_values("n_sig", ascending=False).head(15).to_string(index=False))
print("\n== bottom 10 runs by n_sig ==\n", pr.sort_values("n_sig").head(10).to_string(index=False))
c = t[t.covered].groupby(["event_id", "cell_line"], observed=True).size().unstack(fill_value=0)
nruns = pr.groupby("cell_line").size()
cov = pd.DataFrame({cl: c[cl] / nruns[cl] for cl in nruns.index}).reindex(ev.event_id).fillna(0)
cov = cov.add_prefix("covfrac_").reset_index().merge(ev, on="event_id")
cov.to_csv(f"{U}/summary_event_coverage.tsv", sep="\t", index=False)
print("\n== events covered (both groups >=10 reads) in a fraction of runs, per cell line ==")
for cl in nruns.index:
    for f in (0.5, 0.8, 1.0):
        m = cov[f"covfrac_{cl}"] >= f; print(f"{cl} >= {f:.0%} of {nruns[cl]} runs: {m.sum():,} events", cov[m].type.value_counts().to_dict())
both = (cov.covfrac_K562 >= 0.5) & (cov.covfrac_HepG2 >= 0.5); print("both cell lines >=50%:", both.sum())
s = t[t.covered]; print("\n== covered (event, run) pairs:", len(s), "| |dPSI|>0.1:", round((s.dpsi.abs() > 0.1).mean(), 4), "| sig:", round(s.sig.mean(), 4), "| |dPSI|>0.05:", round((s.dpsi.abs() > 0.05).mean(), 4))
print("dPSI quantiles among covered:", s.dpsi.quantile([.01, .05, .25, .5, .75, .95, .99]).round(3).to_dict())
es = t[t.sig].groupby("event_id", observed=True).size(); print("events significant in >=1 run:", len(es), "| >=5 runs:", (es >= 5).sum(), "| >=20 runs:", (es >= 20).sum(), "| max:", es.max())
