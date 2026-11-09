"""Confirmatory checks on the composition-matched survivors SRSF3 (high reproducibility, SR protein) and CSNK1E (fallopian cohort):
 a matched pairs restricted to the homogeneous stratum (premenopausal, batch C) -> removes batch/menopause as well as composition
 b activity vs each composition score over the whole cohort (is SRSF3 activity itself composition-like?)
 c breast cohort: composition-matched pairs for SRSF3/CSNK1E (replication, sign)
 d event level, regression-free: on SRSF3's own clean events, fraction of events where the carrier-minus-matched-control PSI difference has the
   sign opposite to the SRSF3 knockdown effect (expected > 0.5 if carriers behave like 'more SRSF3'); null = same statistic for 200 random RBPs' clean events
 e independent layer: SRSF3 and CSNK1E mRNA (log2 CPM) carriers vs non-carriers (raw, composition-matched, and full-adjusted)
Output: stdout -> logs/srsf3_followup.txt"""
import sys, numpy as np, pandas as pd; from scipy import stats; sys.path.insert(0, "/well/blagden/users/eic789/PerturbBridge/14_inverse"); import retrieval as R
PB = R.PB; O = R.O; rng = np.random.default_rng(5); pd.set_option("display.width", 260)
compcols = ["epithelial", "ft_secretory", "ft_ciliated", "breast_luminal", "breast_basal", "stroma", "smooth_muscle", "immune", "endothelial", "adipose"]
def matched(Z, grp):
    ci, ni = np.where(grp == 1)[0], np.where(grp == 0)[0]; return [(i, ni[np.argmin(np.sqrt(((Z[ni] - Z[i]) ** 2).sum(1)))]) for i in ci]
def paired(S, pairs, h, n=2000):
    d = np.array([S.iloc[i][h] - S.iloc[j][h] for i, j in pairs]); d = d[np.isfinite(d)]; t, p = stats.ttest_1samp(d, 0); null = [abs(np.mean(d * rng.choice([-1, 1], len(d)))) for _ in range(n)]
    return d.mean(), t, (np.sum(np.array(null) >= abs(d.mean())) + 1) / (n + 1), len(d)
# fallopian
S = pd.read_parquet(f"{O}/h4ext_activity_fallopian_splicing.parquet"); samples = S.index.tolist(); C = f"{PB}/05_precancer_cohorts/GSE274628"
st = pd.read_csv(f"{C}/sample_table.tsv", sep="\t").set_index("gsm").loc[samples]; comp = pd.read_csv(f"{C}/cell_composition_scores.tsv", sep="\t").set_index("gsm").loc[samples]
carrier = (st.genotype != "Non-BRCAm").astype(float).values; Z = ((comp[compcols] - comp[compcols].mean()) / comp[compcols].std()).values
print("== a. matched pairs within premenopausal & batch C ==")
sub = np.where((st.reproductive_status == "Pre").values & (st.batch == "C").values)[0]; Ss = S.iloc[sub]; pr = matched(Z[sub], carrier[sub])
for h in ["SRSF3", "CSNK1E", "RRP9", "EIF3D"]:
    m, t, p, n = paired(Ss, pr, h); print(f"  {h:7s} pairs {n}  mean diff {m:+.4f}  t {t:+.2f}  sign-flip p {p:.3f}")
print("\n== b. activity vs composition (whole fallopian cohort) ==")
for h in ["SRSF3", "CSNK1E"]: print(f"  {h:7s}", {c: round(stats.pearsonr(S[h].values, comp[c].values)[0], 2) for c in compcols}, "| vs carrier r", round(stats.pointbiserialr(carrier, S[h].values)[0], 2))
print("\n== c. breast cohort, composition-matched pairs ==")
Sb = pd.read_parquet(f"{O}/h4ext_activity_breast_splicing.parquet"); Cb = f"{PB}/05_precancer_cohorts/GSE166044"; d = f"{Cb}/rmats_susceptible_vs_healthy"
g1 = [p.split("/")[-1].split(".")[0] for p in open(f"{d}/b1.txt").read().strip().split(",")]; grp_b = np.array([1.0 if s in g1 else 0.0 for s in Sb.index]); compb = pd.read_csv(f"{Cb}/cell_composition_scores.tsv", sep="\t").set_index("gsm").loc[Sb.index]
Zb = ((compb[compcols] - compb[compcols].mean()) / (compb[compcols].std() + 1e-9)).values; prb = matched(Zb, grp_b)
for h in ["SRSF3", "CSNK1E"]:
    m, t, p, n = paired(Sb, prb, h); print(f"  {h:7s} pairs {n}  mean diff {m:+.4f}  t {t:+.2f}  sign-flip p {p:.3f}  (fallopian sign was negative)")
print("\n== d. event level on SRSF3 clean events (fallopian, composition-matched pairs) ==")
T = R.rmats_psi_table(f"{C}/rmats_carrier_vs_noncarrier"); g1f = [p.split("/")[-1].split(".")[0] for p in open(f"{C}/rmats_carrier_vs_noncarrier/b1.txt").read().strip().split(",")]; g2f = [p.split("/")[-1].split(".")[0] for p in open(f"{C}/rmats_carrier_vs_noncarrier/b2.txt").read().strip().split(",")]
T["sample"] = T["sample"].map({f"g1_{i}": s for i, s in enumerate(g1f)} | {f"g2_{i}": s for i, s in enumerate(g2f)}); T.loc[T["cov"] < 10, "psi"] = np.nan
W = T.pivot_table(index="sample", columns="event_id", values="psi", aggfunc="first").reindex(samples); W = W.loc[:, W.notna().mean(0) >= 0.7]
pairs = matched(Z, carrier); Dm = np.nanmean(np.stack([W.iloc[i].values - W.iloc[j].values for i, j in pairs]), 0); Dm = pd.Series(Dm, index=W.columns)
runs = pd.concat([pd.read_parquet(f"{O}/dict_splicing_runs_{cl}.parquet") for cl in ("HepG2", "K562")]); clean = runs[runs.clean].groupby(["target", "event_id"]).dpsi_cc.mean()
def frac_opposite(tg, min_n=30):
    e = clean.loc[tg] if tg in clean.index.get_level_values(0) else None
    if e is None: return np.nan, 0
    e = e[e.index.isin(Dm.index)]; e = e[np.isfinite(Dm[e.index].values)]
    if len(e) < min_n: return np.nan, len(e)
    return float(np.mean(np.sign(Dm[e.index].values) == -np.sign(e.values))), len(e)
f, n = frac_opposite("SRSF3"); print(f"  SRSF3 clean events in cohort: {n}; fraction where carrier-minus-matched-control PSI change is OPPOSITE to SRSF3 knockdown: {f:.3f}")
targets = [t for t in clean.index.get_level_values(0).unique() if t != "SRSF3"]; null = [frac_opposite(t)[0] for t in rng.choice(targets, 200, replace=False)]; null = np.array([x for x in null if np.isfinite(x)])
print(f"  same statistic for {len(null)} random RBPs: mean {null.mean():.3f}, 95th pct {np.percentile(null, 95):.3f}, p = {(np.sum(null >= f) + 1) / (len(null) + 1):.3f}")
w = Dm[clean.loc['SRSF3'].index.intersection(Dm.index)]; e = clean.loc["SRSF3"][w.index]; print(f"  Pearson r(carrier-minus-control dPSI, SRSF3 knockdown dpsi_cc) on these events: {stats.pearsonr(w.values[np.isfinite(w.values)], e.values[np.isfinite(w.values)])[0]:+.3f}")
print("\n== e. SRSF3 / CSNK1E mRNA (log2 CPM) ==")
M = pd.DataFrame({s: pd.read_csv(f"{C}/star/{s}.ReadsPerGene.out.tab", sep="\t", header=None, index_col=0).iloc[4:, 0] for s in samples}).T; M.columns = [c.split(".")[0] for c in M.columns]; L = np.log2(M.div(M.sum(1), axis=0) * 1e6 + 1)
ids = {"SRSF3": "ENSG00000112081", "CSNK1E": "ENSG00000213923"}
pre = (st.reproductive_status == "Pre").astype(float).values; batch = pd.get_dummies(st.batch.astype(str), drop_first=True).values.astype(float)
for g, gid in ids.items():
    if gid not in L.columns: print(f"  {g}: {gid} not in counts"); continue
    y = L[gid].values; a, b = y[carrier == 1], y[carrier == 0]; t, p = stats.ttest_ind(a, b, equal_var=False)
    dm = np.array([y[i] - y[j] for i, j in pairs]); tm, pm = stats.ttest_1samp(dm, 0)
    X = np.hstack([np.ones((len(y), 1)), carrier[:, None], batch, pre[:, None], Z]); bb = np.linalg.lstsq(X, y, rcond=None)[0]; res = y - X @ bb; se = np.sqrt((res ** 2).sum() / (len(y) - X.shape[1]) * np.linalg.pinv(X.T @ X)[1, 1])
    print(f"  {g}: carrier mean {a.mean():.2f} vs {b.mean():.2f} (t {t:+.2f}, p {p:.3f}) | matched-pair diff {dm.mean():+.3f} (p {pm:.3f}) | full-adjusted beta {bb[1]:+.3f} (t {bb[1] / se:+.2f})")
print("DONE")
