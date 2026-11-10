# Programme-level forward model (v5): plan written before any result (2026-11-10)

## Why
Event-level dPSI of one RBP agrees across cell lines with median r 0.05; scores on principal response programmes agree with r 0.52 (programme 1) and 0.21-0.39 (programmes 2-10) for RBPs not used to fit the programmes (13_model/13_lowrank_reproducibility_check.py). The RBP descriptors used so far (mean ESM2, GenePT) barely separate RBPs (median pairwise cosine 0.978) and a ridge model with shuffled RBP identities scores the same as with true identities.

## Design
- Events: covered in >= 80 % of RBPs in both cell lines. Matrix Y[cell]: RBP x event, RBP-mean control-centred dPSI, missing = training-fold event mean.
- Folds: the existing paralog-grouped 5 folds (11_splits/rbp_folds.tsv). Everything that is fitted (event means, programme basis, feature scaling of supervised parts, regression weights, hyper-parameters) uses training-fold RBPs only.
- Target: scores of each RBP on the first k = 10 right singular vectors of the centred training matrix (per cell line and fold).
- Predictor: multi-output ridge from RBP descriptors to scores; lambda by grouped inner cross-validation over the training folds. Secondary predictor: mean score of the 5 nearest training RBPs in descriptor space.
- Descriptor sets, each added to the previous: S0 sequence/text (ESM2 + GenePT, 32 PCs); S1 + genome-wide Perturb-seq expression signature of the same gene in K562 (Replogle et al. 2022; PCs fitted on all perturbations, unsupervised); S2 + STRING v12 network embedding and Complex Portal membership; S3 + DepMap 24Q4 CRISPR gene-effect profile (PCs across genes); S4 + observational prior: correlation of RBP expression with event PSI across TCGA tumours (within project, Fisher-averaged), projected on the training programme basis.
- Controls for every set: descriptors permuted across RBPs (20 permutations); zero predictor (training mean); cross-cell-line oracle (the same RBP's measured scores in the other cell line), which is the reference ceiling.

## Metrics (held-out RBPs, pooled over the 5 folds)
1. Primary: programme-1 score Pearson r, and subspace R2 over the 10 programmes (1 - SSE/SST against the training mean).
2. Per-programme Pearson r for programmes 2-10 (fold bases aligned to a reference basis by absolute cosine).
3. Event level, for comparison with the v3/ridge benchmark: per-RBP Pearson r between reconstructed and measured centred profiles; direction accuracy on clean (RBP, event) pairs for the reconstruction, against the training-mean prior alone.
4. Stratified by reproducibility class (14_inverse/dict_rbps.tsv).

## Decision rules (fixed now)
- A descriptor set "works" if its primary metrics exceed the 95th percentile of its 20 permuted controls in both cell lines, and the paired bootstrap 95 % CI (2,000 resamples over RBPs) of model minus permuted control excludes 0.
- "Qualitative change" is claimed only if programme-1 r reaches at least half of the cross-cell-line oracle r in both cell lines AND at least 3 of programmes 2-10 pass the rule above.
- An added set is kept only if it improves subspace R2 over the previous set with a bootstrap CI excluding 0 in at least one cell line and no significant loss in the other.
- If S1-S4 fail, the result is reported as negative; no further descriptor search on these folds.

## Scope note
With S1 the claim changes from "predict from sequence alone" to "predict the splicing response of an RBP never profiled for splicing from its expression-layer perturbation signature and network context". Folds stay sequence-paralog grouped; functional neighbours (same complex) may sit in training folds, which is the intended source of transfer and is reported as such.
