"""Block 4a (post-freeze in-house validation): HeyA8 LARP1 knockout vs wild type rMATS (nf-core rnasplice run,
LncRNASPLICE/output_lnc/star_salmon/rmats/KO-WT/rmats_post), mapped onto the ENCODE dictionary. LARP1 itself is not an ENCODE
knockdown target; LARP4 and LARP7 are the nearest family members and are the pre-registered expectation (DESIGN_v2 week 14).
Query = events with nominal p < 0.01 and |dPSI| > 0.1 (coverage >= 10 both groups), capped at 500 by p (3 vs 3 replicates give only 66 FDR-significant events); min shared events 30. Score and rank as in H3/H4;
null = 200 random event subsets of the same size drawn from all tested events with their observed dPSI (what the dictionary scores
for an unselected profile). Outputs: 14_inverse/heya8_ranking.tsv; stdout -> logs/heya8_validation.txt"""
import sys, numpy as np, pandas as pd; sys.path.insert(0, "/well/blagden/users/eic789/PerturbBridge/14_inverse"); import retrieval as R
PB = R.PB; O = R.O; D = "/well/blagden/users/eic789/LncRNASPLICE/output_lnc/star_salmon/rmats/KO-WT/rmats_post"; rng = np.random.default_rng(0); pd.set_option("display.width", 250)
P = R.load_dictionary(); info = R.rbp_info(); ev = pd.read_parquet(f"{R.U}/events.parquet", columns=["key", "event_id"]); k2i = dict(zip(ev.key, ev.event_id)); rows = []
for t, cols in R.COORD.items():
    x = pd.read_csv(f"{D}/{t}.MATS.JC.txt", sep="\t", dtype=str); key = t + "|" + x["chr"] + "|" + x["strand"] + "|" + x[cols].astype(str).agg("|".join, axis=1); x["event_id"] = key.map(k2i); x = x[x.event_id.notna()]
    c1 = x.IJC_SAMPLE_1.str.split(",", expand=True).astype(float).values + x.SJC_SAMPLE_1.str.split(",", expand=True).astype(float).values; c2 = x.IJC_SAMPLE_2.str.split(",", expand=True).astype(float).values + x.SJC_SAMPLE_2.str.split(",", expand=True).astype(float).values
    rows.append(pd.DataFrame({"event_id": x.event_id.values, "type": t, "dpsi": x.IncLevelDifference.astype(float).values, "fdr": x.FDR.astype(float).values, "p": x.PValue.astype(float).values, "cov1": c1.mean(1), "cov2": c2.mean(1), "n1": c1.shape[1], "n2": c2.shape[1]}))
T = pd.concat(rows); T = T[(T.cov1 >= 10) & (T.cov2 >= 10) & np.isfinite(T.dpsi)].drop_duplicates("event_id").set_index("event_id")
print("HeyA8 KO vs WT: group sizes", T.n1.iloc[0], T.n2.iloc[0], "| tested events in universe", len(T), "| FDR<0.05 & |dPSI|>0.1:", int(((T.fdr < 0.05) & (T.dpsi.abs() > 0.1)).sum()), "| sign convention: rMATS b1 - b2 (check bamlist: KO first?)")
sig = T[(T.p < 0.01) & (T.dpsi.abs() > 0.1)].sort_values("p").head(500); q = sig.dpsi; print("query: p<0.01 & |dPSI|>0.1 ->", len(q), "events (FDR<0.05 gave only 66)")
s, n = R.score_against(q, P, min_shared=30); rk = R.rank_report(s, info, "HeyA8_LARP1_KO")
null = np.full((200, len(P.columns)), np.nan)
for k in range(200):
    sub = T.sample(len(q), random_state=int(rng.integers(1e9))).dpsi; sk, _ = R.score_against(sub, P, min_shared=30); null[k] = sk.reindex(P.columns).values
obs = s.reindex(P.columns).values; rk = rk.merge(pd.DataFrame({"rbp": P.columns, "perm_p": [(np.nansum(null[:, j] >= obs[j]) + 1) / (np.isfinite(null[:, j]).sum() + 1) if np.isfinite(obs[j]) else np.nan for j in range(len(obs))], "null_mean": np.nanmean(null, 0), "null_sd": np.nanstd(null, 0)}), on="rbp")
rk["z_vs_null"] = (rk.score - rk.null_mean) / rk.null_sd; rk.to_csv(f"{O}/heya8_ranking.tsv", sep="\t", index=False, float_format="%.4f")
print(f"query events {len(q)}, RBPs scored {s.notna().sum()}\n"); print(rk.head(20).round(3).to_string(index=False))
print("\n== LARP family and other pre-registered comparators =="); print(rk[rk.rbp.isin(["LARP4", "LARP7", "LARP1", "PABPC4", "PABPC1", "HNRNPA1"])].round(3).to_string(index=False))
print("\nbamlist:", open("/well/blagden/users/eic789/LncRNASPLICE/output_lnc/star_salmon/rmats/bamlist").read()[:400] if __import__("os").path.isfile("/well/blagden/users/eic789/LncRNASPLICE/output_lnc/star_salmon/rmats/bamlist") else __import__("os").listdir("/well/blagden/users/eic789/LncRNASPLICE/output_lnc/star_salmon/rmats/bamlist"))
print("DONE", f"{O}/heya8_ranking.tsv")
