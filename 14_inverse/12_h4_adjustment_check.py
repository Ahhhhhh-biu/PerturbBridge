"""Is the full composition adjustment an over-correction? Four direct checks on the fallopian cohort (activity matrices from 11_h4_robustness.py):
 1 Do carriers and non-carriers differ in cell composition? (per-score t test; R2 of carrier status explained by the 10 scores)
 2 Does the RBP activity track composition WITHIN non-carriers only (no carrier status involved)? If yes, the activity measures composition.
 3 Graded adjustment: none -> batch+menopause -> +epithelial -> +ciliated,stroma,smooth muscle,endothelial -> all 10 (+age). Family-wise p and t of the
   five second-stage hits at each level, so the reader sees where the signal goes.
 4 Regression-free check: each carrier matched to its nearest non-carrier in composition space (10 standardised scores); paired activity difference,
   sign-flip permutation. No covariates are 'regressed out' here, so no over-correction is possible.
Output: 14_inverse/h4adj_check.tsv (graded table); stdout -> logs/h4_adjustment_check.txt"""
import sys, numpy as np, pandas as pd; from scipy import stats; sys.path.insert(0, "/well/blagden/users/eic789/PerturbBridge/14_inverse"); import retrieval as R
PB = R.PB; O = R.O; rng = np.random.default_rng(3); pd.set_option("display.width", 260); NPERM = 1000; HITS = ["RRP9", "EIF3D", "EIF2S1", "SRSF3", "SUGP2", "TAF15"]
S = pd.read_parquet(f"{O}/h4ext_activity_fallopian_splicing.parquet"); samples = S.index.tolist()
C = f"{PB}/05_precancer_cohorts/GSE274628"; st = pd.read_csv(f"{C}/sample_table.tsv", sep="\t").set_index("gsm").loc[samples]; comp = pd.read_csv(f"{C}/cell_composition_scores.tsv", sep="\t").set_index("gsm").loc[samples]
compcols = ["epithelial", "ft_secretory", "ft_ciliated", "breast_luminal", "breast_basal", "stroma", "smooth_muscle", "immune", "endothelial", "adipose"]
carrier = (st.genotype != "Non-BRCAm").astype(float).values; pre = (st.reproductive_status == "Pre").astype(float).values; batch = pd.get_dummies(st.batch.astype(str), drop_first=True).values.astype(float); age = pd.to_numeric(st.age_at_time_of_surgery, errors="coerce"); age = age.fillna(age.median()).values
Z = (comp[compcols] - comp[compcols].mean()) / comp[compcols].std()
print("== 1. composition differs between carriers and non-carriers? ==")
for c in compcols:
    a, b = comp.loc[carrier == 1, c], comp.loc[carrier == 0, c]; t, p = stats.ttest_ind(a, b, equal_var=False); print(f"  {c:15s} carrier mean {a.mean():6.2f}  non-carrier {b.mean():6.2f}  t {t:+5.2f}  p {p:.3g}")
X = np.hstack([np.ones((len(carrier), 1)), Z.values]); b_ = np.linalg.lstsq(X, carrier, rcond=None)[0]; r2 = 1 - ((carrier - X @ b_) ** 2).sum() / ((carrier - carrier.mean()) ** 2).sum()
Xb = np.hstack([np.ones((len(carrier), 1)), batch, pre[:, None]]); b2 = np.linalg.lstsq(Xb, carrier, rcond=None)[0]; r2b = 1 - ((carrier - Xb @ b2) ** 2).sum() / ((carrier - carrier.mean()) ** 2).sum()
print(f"  R2 of carrier status explained by the 10 composition scores: {r2:.2f} | by batch + menopause: {r2b:.2f}")
print("\n== 2. within NON-CARRIERS only: correlation of hit-RBP activity with composition ==")
nc = carrier == 0
for h in HITS:
    v = S.loc[nc, h].values; rr = {c: round(stats.pearsonr(v, comp.loc[nc, c].values)[0], 2) for c in ["ft_ciliated", "ft_secretory", "epithelial", "stroma", "smooth_muscle", "endothelial", "immune"]}; print(f"  {h:7s} {rr}")
print("\n== 3. graded adjustment (family-wise p of max |t| over 478 RBPs; t of the hits) ==")
levels = {"none": [], "batch+menopause": [batch, pre[:, None]], "+epithelial": [batch, pre[:, None], Z[["epithelial"]].values], "+ciliated,stroma,smooth muscle,endothelial": [batch, pre[:, None], Z[["epithelial", "ft_ciliated", "stroma", "smooth_muscle", "endothelial"]].values], "all 10 + age": [batch, pre[:, None], Z.values, ((age - age.mean()) / age.std())[:, None]]}
rows = []
for name, cov in levels.items():
    X = np.hstack([np.ones((len(carrier), 1)), carrier[:, None]] + cov); b, t, p = R.ols_group_effect(S.values, X); tm = np.nanmax(np.abs(t))
    null = np.array([np.nanmax(np.abs(R.ols_group_effect(S.values, np.hstack([np.ones((len(carrier), 1)), rng.permutation(carrier)[:, None]] + cov))[1])) for _ in range(300)])
    pf = (np.sum(null >= tm) + 1) / 301; th = {h: round(float(t[list(S.columns).index(h)]), 2) for h in HITS}; rows.append(dict(adjustment=name, n_covariates=X.shape[1] - 2, familywise_p=pf, max_abs_t=round(tm, 2), **{f"t_{h}": v for h, v in th.items()}))
    print(f"  {name:45s} cov {X.shape[1] - 2:2d} | family-wise p {pf:.3f} | max|t| {tm:.2f} | {th}")
pd.DataFrame(rows).to_csv(f"{O}/h4adj_check.tsv", sep="\t", index=False)
print("\n== 4. composition-matched pairs (no regression) ==")
ci, ni = np.where(carrier == 1)[0], np.where(carrier == 0)[0]; Zv = Z.values; pairs = []
for i in ci:
    d = np.sqrt(((Zv[ni] - Zv[i]) ** 2).sum(1)); j = ni[np.argmin(d)]; pairs.append((i, j, d.min()))
print(f"  {len(pairs)} carriers each matched to nearest non-carrier (median composition distance {np.median([p[2] for p in pairs]):.2f} SD units; unmatched carrier-vs-all median {np.median([np.sqrt(((Zv[ni] - Zv[i]) ** 2).sum(1)).mean() for i in ci]):.2f})")
for h in HITS + ["all 478 RBPs: max |mean paired diff| / sd"]:
    if h.startswith("all"):
        D = np.array([S.iloc[i].values - S.iloc[j].values for i, j, _ in pairs]); tt = np.nanmean(D, 0) / (np.nanstd(D, 0, ddof=1) / np.sqrt(len(D))); tm = np.nanmax(np.abs(tt))
        null = [np.nanmax(np.abs(np.nanmean(D * rng.choice([-1, 1], (len(D), 1)), 0) / (np.nanstd(D, 0, ddof=1) / np.sqrt(len(D))))) for _ in range(NPERM)]; print(f"  family-wise: max paired |t| {tm:.2f}, sign-flip p {(np.sum(np.array(null) >= tm) + 1) / (NPERM + 1):.3f}, top RBP {S.columns[np.nanargmax(np.abs(tt))]}"); continue
    d = np.array([S.iloc[i][h] - S.iloc[j][h] for i, j, _ in pairs]); d = d[np.isfinite(d)]; t, p = stats.ttest_1samp(d, 0)
    null = [abs(np.mean(d * rng.choice([-1, 1], len(d)))) for _ in range(NPERM)]; print(f"  {h:7s} mean paired diff {d.mean():+.4f}  t {t:+.2f}  parametric p {p:.3f}  sign-flip p {(np.sum(np.array(null) >= abs(d.mean())) + 1) / (NPERM + 1):.3f}")
print("DONE")
