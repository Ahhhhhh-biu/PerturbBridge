"""Aggregate held-out results of the forward model over folds and seeds (13_model/runs/<cell>_f<k>_<tier>[_s<seed>]/per_run.tsv),
compare with the control-centred baselines (12_baselines/baselines_<cell>_dpsi_cc_per_run.tsv) and stratify by data tier, by the
held-out RBP's number of significant events (>=1000 / 100-1000 / <100) and by its closest training paralog (BLASTP max identity with
coverage >= 0.3 to any RBP in another fold: 35-40 / 30-35 / <30; the paralog-grouped split puts every pair >= 40 % in one fold).
Unit = RBP (mean over its runs and seeds); CI = bootstrap over RBPs. Both the plain RBP mean and the n_sig-weighted mean are
reported, because most RBPs have few significant events and their per-RBP metrics are noise.
Outputs 13_model/results_summary.tsv, results_by_rbp.tsv, results_by_seed.tsv (+ stdout)."""
import glob, os, re, numpy as np, pandas as pd
TAG = os.environ.get("TAG", "")  # e.g. TAG=_v2 aggregates runs/<cell>_f<k>_<tier>_v2_s<seed>/ and writes results_*_v2.tsv
PB = "/well/blagden/users/eic789/PerturbBridge"; rng = np.random.default_rng(0); met = ["r_all", "r_sig", "sign_sig", "sign_0.05_0.1", "sign_0.1_0.2", "sign_0.2_1.01", "mse_all"]
folds = pd.read_csv(f"{PB}/11_splits/rbp_folds.tsv", sep="\t").set_index("target"); bl = pd.read_csv(f"{PB}/11_splits/blast/allvsall.tsv", sep="\t", names=["q", "s", "pident", "length", "qlen", "slen", "evalue", "bits"])
bl = bl[(bl.q != bl.s) & (bl.length >= 0.3 * np.minimum(bl.qlen, bl.slen))]; bl["qf"] = bl.q.map(folds.fold); bl["sf"] = bl.s.map(folds.fold)
near = bl[bl.qf != bl.sf].groupby("q").pident.max(); strat = lambda t: "35-40" if near.get(t, 0) >= 35 else ("30-35" if near.get(t, 0) >= 30 else "<30")
R = []
files = glob.glob(f"{PB}/13_model/runs/*_f[0-4]_[ABC]{TAG}_s[0-9]/per_run.tsv") if TAG else glob.glob(f"{PB}/13_model/runs/*_f[0-4]_[ABC]/per_run.tsv") + glob.glob(f"{PB}/13_model/runs/*_f[0-4]_[ABC]_s[0-9]/per_run.tsv")
print("run files:", len(files))
for f in files:
    d = pd.read_csv(f, sep="\t"); dn = f.split("/")[-2]; d["cell_line"] = dn.split("_")[0]; d["method"] = "model_" + d.tier_model
    m = re.search(r"_s(\d+)$", dn); d["seed"] = int(m.group(1)) if m else 0; R.append(d)
for cl in ("HepG2", "K562"):
    b = pd.read_csv(f"{PB}/12_baselines/baselines_{cl}_dpsi_cc_per_run.tsv", sep="\t"); b["cell_line"] = cl; b["tier_data"] = b.tier; b["method"] = "baseline_" + b.method; b["seed"] = 0; R.append(b)
R = pd.concat(R, ignore_index=True); R["paralog_strat"] = R.target.map(strat)
nseed = R.groupby("method").seed.nunique(); print("seeds per method:", nseed.to_dict())
# per RBP: mean over runs and seeds; n_sig = significant events summed over runs (not over seeds)
by = R.groupby(["cell_line", "method", "target"]).agg({**{m: "mean" for m in met}, "tier_data": "first", "paralog_strat": "first"}).reset_index()
ns = R.drop_duplicates(["cell_line", "method", "target", "run"]).groupby(["cell_line", "method", "target"]).n_sig.sum().rename("n_sig").reset_index()
by = by.merge(ns, on=["cell_line", "method", "target"]); by["nsig_strat"] = pd.cut(by.n_sig, [-1, 99, 999, 10**9], labels=["<100", "100-1000", ">=1000"]).astype(str)
by.to_csv(f"{PB}/13_model/results_by_rbp{TAG}.tsv", sep="\t", index=False)
# per seed (model only): RBP-mean and n_sig-weighted sign_sig / r_sig, to show seed-to-seed spread
bs_ = R[R.method.str.startswith("model_")].groupby(["cell_line", "method", "seed", "target"]).agg(r_sig=("r_sig", "mean"), sign_sig=("sign_sig", "mean"), n_sig=("n_sig", "sum")).reset_index()
seed_tab = bs_.groupby(["cell_line", "method", "seed"]).apply(lambda g: pd.Series({"r_sig": g.r_sig.mean(), "sign_sig": g.sign_sig.mean(), "r_sig_w": np.average(g.r_sig.fillna(0), weights=g.n_sig), "sign_sig_w": np.average(g.sign_sig.fillna(0), weights=g.n_sig)})).reset_index()
seed_tab.to_csv(f"{PB}/13_model/results_by_seed{TAG}.tsv", sep="\t", index=False)
def boot(g):
    out = {"n_rbp": len(g)}
    for m in met:
        v = g[m].dropna().values
        if len(v) < 3: out[m] = "NA"; continue
        bs = [rng.choice(v, len(v)).mean() for _ in range(2000)]; out[m] = f"{v.mean():.3f} [{np.percentile(bs, 2.5):.3f},{np.percentile(bs, 97.5):.3f}]"
    for m in ("r_sig", "sign_sig"):  # n_sig-weighted mean over RBPs, bootstrap over RBPs
        gg = g.dropna(subset=[m]); v = gg[m].values; w = gg.n_sig.values.astype(float)
        if len(v) < 3 or w.sum() == 0: out[m + "_w"] = "NA"; continue
        bs = []
        for _ in range(2000):
            i = rng.integers(0, len(v), len(v)); bs.append(np.average(v[i], weights=w[i]))
        out[m + "_w"] = f"{np.average(v, weights=w):.3f} [{np.percentile(bs, 2.5):.3f},{np.percentile(bs, 97.5):.3f}]"
    return pd.Series(out)
S = []
for name, keys in [("all", ["cell_line", "method"]), ("by_data_tier", ["cell_line", "method", "tier_data"]), ("by_nsig", ["cell_line", "method", "nsig_strat"]), ("by_paralog", ["cell_line", "method", "paralog_strat"])]:
    s = by.groupby(keys).apply(boot).reset_index(); s.insert(0, "stratum", name); S.append(s)
S = pd.concat(S, ignore_index=True); S.to_csv(f"{PB}/13_model/results_summary{TAG}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
print("\n== per seed (model only) ==\n", seed_tab.round(3).to_string(index=False))
for st in ("all", "by_data_tier", "by_nsig", "by_paralog"): print(f"\n== {st} ==\n", S[S.stratum == st].drop(columns=["stratum", "sign_0.05_0.1", "sign_0.1_0.2", "r_all"]).to_string(index=False))
# paired comparison model A vs B vs C on the same RBPs (H1), and vs best baseline; plain and n_sig-weighted
def paired(sub, a_, b_, m, weighted):
    w = sub.pivot(index="target", columns="method", values=m); n = sub.pivot(index="target", columns="method", values="n_sig")
    if a_ not in w or b_ not in w: return None
    d = (w[a_] - w[b_]).dropna(); wt = n[a_].reindex(d.index).fillna(0).values.astype(float) if weighted else np.ones(len(d))
    bs = []
    for _ in range(2000):
        i = rng.integers(0, len(d), len(d)); bs.append(np.average(d.values[i], weights=wt[i]))
    return f"{np.average(d.values, weights=wt):+.4f} [{np.percentile(bs, 2.5):+.4f},{np.percentile(bs, 97.5):+.4f}] (n={len(d)})"
print("\n== paired differences (sign_sig | r_sig), RBP mean and n_sig-weighted ==")
for cl in ("HepG2", "K562"):
    sub = by[by.cell_line == cl]
    for a_, b_ in [("model_A", "model_B"), ("model_B", "model_C"), ("model_A", "model_C"), ("model_C", "baseline_ridge")]:
        print(f"{cl} {a_} - {b_}: sign {paired(sub, a_, b_, 'sign_sig', False)} | weighted {paired(sub, a_, b_, 'sign_sig', True)} | r_sig {paired(sub, a_, b_, 'r_sig', False)} | weighted {paired(sub, a_, b_, 'r_sig', True)}")
    tA = sub[sub.tier_data == "A"]
    print(f"{cl} H1 core (tier-A RBPs) model_A - model_B: sign {paired(tA, 'model_A', 'model_B', 'sign_sig', False)} | weighted {paired(tA, 'model_A', 'model_B', 'sign_sig', True)} | r_sig {paired(tA, 'model_A', 'model_B', 'r_sig', False)} | weighted {paired(tA, 'model_A', 'model_B', 'r_sig', True)}")
    big = sub[sub.n_sig >= 1000]
    print(f"{cl} RBPs with >=1000 sig events model_A - model_B: sign {paired(big, 'model_A', 'model_B', 'sign_sig', False)} | model_A - model_C: {paired(big, 'model_A', 'model_C', 'sign_sig', False)}")
