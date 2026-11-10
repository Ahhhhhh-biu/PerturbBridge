# Second batch of cohorts for the inverse mapping: plan written before any result (2026-11-10)

Cohorts were chosen from GEO by design criteria learnt from the first four cohorts (raw reads public, fresh or sorted material, >= 20 M reads, >= 20 per group where possible). Rejected after inspection: GSE281303 (DCIS progressors; no raw data), GSE232988 (myeloma precursor states; no raw data, single-cell), GSE292971 (Barrett's, single-cell), GSE302239 (oral, single-end, median 9.8 M reads), GSE211032 (monosomy 7; not a splicing-factor design).

Processing is identical to the earlier cohorts: fastp, STAR (GRCh38, GENCODE v47), rMATS-turbo, events keyed by coordinates to the ENCODE universe; per-event OLS of PSI on group and covariates (coverage >= 10, events measured in >= 70 % of samples); query = p < 0.05 and |beta| > 0.05, at most 500 events; Pearson r against the 478-RBP dictionary; 200 label permutations repeating the whole procedure; per-RBP permutation p, BH over RBPs, and p_top = P(maximum of the null >= maximum observed).

## A. Known-answer test in a pre-leukaemic clonal disease: GSE220518 (bone-marrow CD34+ cells, lower-risk MDS; 185 patients, 7 healthy)
Somatic splicing-factor mutations give tests with a known regulator in patient cells that are already purified.
- A1 SF3B1-mutant (only SF3B1 mutated) versus MDS without SF3B1/SRSF2/U2AF1 mutation. Expected: the core-spliceosome class. Success = SF3B1 in the top 10 % of the dictionary, or at least three SF3A/SF3B members in the top ten, with p_top < 0.05.
- A2 U2AF1-mutant versus the same wild-type group (U2AF1 is in the dictionary; 9-11 carriers, low power). Success = U2AF1 or U2AF2 in the top 10 % with permutation p < 0.05.
- A3 SRSF2-mutant versus wild type. SRSF2 is not in the dictionary; reported descriptively (which SR proteins, if any, are retrieved).
- Each is run unadjusted and adjusted for five haematopoietic lineage marker scores (stem, erythroid, myeloid, lymphoid, megakaryocytic), because SF3B1-mutant MDS is erythroid-skewed.
## B. Progression question in the same cohort
- B1 patients who progressed to leukaemia (30) versus those who did not (139), adjusted for the three mutation indicators, and additionally for lineage scores. B2 the same within patients without splicing-factor mutations. Negative-by-default criteria as in H4: p_top < 0.05 and more RBPs at permutation p < 0.05 than expected.
- Descriptive: MDS versus healthy CD34+ (7 controls).
## C. Epithelial premalignant series (H4-type, same criteria as before)
- GSE193946 Barrett's oesophagus: high-grade dysplasia (21) vs metaplasia/low-grade (25); adenocarcinoma (11) vs low-grade; ordinal trend. Unadjusted and adjusted for epithelial, stromal and immune marker scores.
- GSE79209 bronchial brushings of smokers: dysplasia (50) vs normal histology (25), adjusted for smoking status, sex and COPD; plus epithelial/ciliated/basal/immune marker scores; ordinal trend over the worst histology grade.
- GSE227919 oral mucosa: premalignant lesions (39) vs control (18); carcinoma (9) vs control; ordinal trend. Unadjusted and composition-adjusted.
- GSE202048 tongue: lesion versus matched margin within patient (OPMD 9 pairs, carcinoma 11 pairs, all 20 pairs); paired differences with sign-flip permutations; unadjusted and adjusted for the paired difference in marker scores.
## Reading rules
- A is a positive control for the method in patient material; C and B are the premalignant questions. A positive A with negative B/C localises the limitation to effect size in tissue rather than to the method.
- Composition-adjusted results take precedence when they disagree with unadjusted ones; a call must survive adjustment to count.
- No contrast definitions, thresholds or covariate sets are changed after the results are seen; additional analyses are labelled post hoc.
