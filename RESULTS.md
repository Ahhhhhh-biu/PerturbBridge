# PerturbBridge: results on one page

Status 2026-11-09. Every number below comes from a file in this repository; the file is named in brackets.

## Question

Given the ENCODE knockdown compendium (HepG2 and K562, 876 knockdown/control pairs, 472 RNA-binding proteins), can a model predict the splicing response (ΔPSI per event) to knocking down an RBP it has never seen? And can the measured responses be used in reverse, to name the perturbed RBP from a splicing shift observed in tissue?

## Setup

- Target: control-centred ΔPSI on a coordinate-keyed event universe (546,917 events; rMATS junction counts). Events whose control PSI deviates from the cell-line consensus are flagged ("control-deviation rule") and excluded from scoring ("clean rows").
- Official split: paralog-grouped 5-fold hold-out of RBPs (BLASTP ≥ 40 % identity components stay in one fold); every held-out RBP has no relative above 40 % identity in its training folds [11_splits/].
- Repeat check: same-RBP independent experiments agree on direction 0.51–0.80 (median 0.56); two thirds of knockdowns have no reproducible splicing effect [12_baselines/].
- Models, all on the same splits and the same clean rows: ridge regression (closed form, event + RBP + interaction blocks, λ = 100), v3 deep model (event encoder × RBP encoder with low-rank bilinear interaction), v4 1D-CNN (900-nt sequence encoder in place of the v3 event encoder). Each model has a shuffled-RBP control trained with permuted RBP identities; v3 also has an event-only control [13_model/benchmark/].

## Results table (RBP-mean on clean held-out rows; HepG2 / K562)

| Model | Direction accuracy | Spearman | MAE clean | MAE on null rows |
|---|---|---|---|---|
| v3 deep | 0.514 / 0.520 | 0.029 / 0.036 | 0.156 / 0.154 | 0.044 / 0.047 |
| v3 shuffled RBP | 0.512 / 0.500 | 0.027 / −0.003 | 0.156 / 0.159 | 0.049 / 0.049 |
| v3 event-only | 0.503 / 0.505 | 0.037 / 0.030 | 0.157 / 0.156 | 0.030 / 0.036 |
| ridge | 0.518 / 0.508 | 0.100 / 0.072 | 0.155 / 0.155 | 0.036 / 0.039 |
| ridge shuffled RBP | 0.517 / 0.508 | 0.100 / 0.076 | 0.156 / 0.156 | 0.037 / 0.040 |
| v4 CNN | 0.516 / 0.523 | 0.034 / 0.045 | 0.155 / 0.153 | 0.041 / 0.042 |
| v4 CNN shuffled RBP | 0.513 / 0.501 | 0.041 / −0.001 | 0.155 / 0.158 | 0.043 / 0.044 |

Paired RBP-level differences with bootstrap 95 % CI [paired_differences.tsv]: v3 − ridge direction −0.004 [−0.013, +0.005] / +0.012 [+0.002, +0.023], Spearman −0.071 [−0.089, −0.053] / −0.036 [−0.059, −0.014]. CNN − v3: +0.002 / +0.003, both n.s. CNN − CNN-shuffled: +0.003 [−0.007, +0.012] / +0.022 [+0.009, +0.035]. Ridge − ridge-shuffled: 0.000 / 0.000.

Ceilings for context [12_baselines/, 13_model/]: same-RBP independent experiments agree on direction 0.56 (median); the zero-parameter control-deviation rule alone reaches 0.60–0.70 on "significant" events when control PSI is given as input, which is why control PSI is not an input here.

## Error analysis (fig_benchmark.png)

- Which RBPs are hard: accuracy tracks the RBP's own reproducibility. In the > 0.7 bin (31 RBPs, mostly core spliceosome) all models reach 0.62–0.72; below 0.6 all are at 0.49–0.52. In that top bin ridge-shuffled (0.673 / 0.721) and event-only (0.557 / 0.685) are as good as v3 (0.678 / 0.694), so most of the accuracy there is an RBP-independent event prior. The RBP-specific gain over the event-only control is +0.12 in HepG2 (n = 11) and ≈ 0 in K562.
- Distance to the training set: max BLASTP identity to any training-fold RBP (20–30 / 30–40 / no hit) does not predict accuracy; all models are flat at 0.49–0.53.
- Calibration: on null rows (events the knockdown did not change) MAE is 0.03–0.05; the deep models predict amplitudes about twice as large as ridge and are the worse calibrated.
- Sequence-specific RBPs with reproducible knockdowns (PTBP1, SRSF1, HNRNPC) are not predicted better than their shuffled controls; eCLIP and motif features add < 0.7 points (tiers A / B / C).

## Programme-level model (v5) [13_model/programme/]

Plan and decision rules were written before the results (PLAN.md). Target: scores of each RBP on 10 response programmes (principal axes of the centred RBP × event matrix, fitted on training-fold RBPs only; 21,077 events covered in both cell lines). Predictor: ridge from RBP descriptors; controls: descriptors permuted across RBPs (20 times) and a cross-cell-line oracle (the same RBP's measured scores in the other cell line).

| Held-out RBPs, HepG2 / K562 | Programme-1 r | Spearman | r without top 5 % | Subspace R² | Programmes 2–10 above control |
|---|---|---|---|---|---|
| permuted descriptors (95th pct) | 0.07 / 0.08 | 0.06 / 0.07 | 0.03 / 0.06 | 0.000 / 0.002 | – |
| sequence + text (as in v3) | 0.35 / 0.43 | 0.09 / 0.15 | 0.03 / 0.15 | 0.023 / 0.076 | 2 / 5 |
| + network (selected by the rule) | 0.45 / 0.60 | 0.19 / 0.20 | 0.14 / 0.26 | 0.060 / 0.176 | 5 / 9 |
| + network + DepMap | 0.51 / 0.60 | 0.25 / 0.24 | 0.26 / 0.31 | 0.074 / 0.169 | 2 / 9 |
| oracle: measured in the other cell line | 0.54 / 0.54 | 0.30 / 0.29 | 0.30 / 0.26 | 0.098 / 0.135 | – |

- Changing the target makes the transferable signal measurable: the same sequence/text descriptors that were indistinguishable from shuffled at the event level predict the programme-1 score with r 0.35 / 0.43 (model minus permuted +0.36 [0.18, 0.51] / +0.43 [0.29, 0.55]).
- Descriptor blocks, added cumulatively (paired bootstrap of subspace R²): K562 Perturb-seq signature −0.010 [−0.022, +0.002] / −0.009 [−0.036, +0.017], not kept; STRING network + Complex Portal +0.028 [+0.009, +0.048] / +0.101 [+0.038, +0.156], kept; DepMap +0.008 [−0.003, +0.018] / +0.010 [−0.022, +0.046], not kept; TCGA expression × PSI co-variation −0.003 / −0.002, not kept (alone: R² −0.003 / 0.001).
- With the selected set the model reaches the oracle in K562 and about half of it in HepG2. Both pre-specified conditions for a qualitative change are met (programme-1 r at least half of the oracle; at least three further programmes above control in both cell lines).
- Caveats. The Pearson r of programme 1 is driven by nine core-spliceosome RBPs with |z| > 3; rank and trimmed correlations are 0.14–0.26. Programme 1 is the core-spliceosome axis (U2AF1/2, SF3A/SF3B, AQR, CDC5L; retained introns 29–40 % of its top events against 12 % overall; loadings agree between cell lines, |cosine| 0.78); programme 2 is shared (0.60); programmes 3–10 are cell-line specific and dominated by single strong knockdowns. Reconstructed event-level direction accuracy is 0.50 / 0.53 (permuted 0.44; ceiling with the 10 measured scores 0.70 / 0.73), so single-event prediction remains weak.

Figure: 13_model/programme/fig_programme.png. Tables: results_summary_final.tsv, results_programme_bootstrap_final.tsv, results_strata_final.tsv, programme_annotation.txt.

## Inverse mapping [14_inverse/]

- H3, ENCODE retrieval across cell lines (query HepG2 profile against the K562 dictionary and vice versa, 406 / 444 candidates): top-5 0.255 / 0.249 versus 0.011 expected; for the high-reproducibility class median rank 1, top-1 0.82 / 0.56. Same-technology versus cross-technology (shRNA vs CRISPR) top-5 0.334 / 0.260 vs 0.237 / 0.256, so most of the signal is biology, not reagent [h3_encode_summary.tsv, h3_cross_technology.tsv].
- H3, TCGA tumours (Snaptron junction counts, cBioPortal mutations; 52 mutant-vs-wild-type tests in five projects, 45 with the gene in the dictionary): negative overall, median rank fraction 0.44 against 0.50 under a uniform null, 4 of 45 in the top 10 % (expected 4.5). The exception is the SF3B1 hotspot in uveal melanoma (17 mutant tumours), whose splicing shift retrieves SF3B3, SF3B1, RBM25 and SF3A3 at the top, the same core-spliceosome class that is nameable in ENCODE. SF3B1 and U2AF1 hotspots in breast, melanoma, lung and leukaemia are not retrieved [h3_tcga_tests.tsv].
- H4, premalignant cohorts (BRCA-carrier fallopian tube GSE274628; pre-diagnosis breast GSE166044): pre-registered group-level test negative (p_top 0.29 / 0.12 against the maximum of the null). A second-stage activity-score analysis produced calls that disappeared under cell-composition adjustment; the composition-matched residue (SRSF3, CSNK1E) is weak, opposite in sign between the two cohorts and without mRNA support. Status: suggestive, unreplicated. A third cohort of FACS-sorted mammary epithelial cells (GSE182338; 120 samples, 48 strains, 24 carrier samples; strain-level permutations) is negative after batch adjustment in every contrast (p_top 0.09–0.91, no BH q < 0.1); its one structured residue is the core-spliceosome block (SF3B1, SF3B3, SF3A1, SF3A3, SF3A2, WBP11, CDC5L) at the top of the adjusted myoepithelial carrier-vs-control ranking (p_top 0.26). A colorectal adenoma → carcinoma cohort (GSE289499; FFPE, Ion Torrent, < 10 M reads) is not evaluable at the splicing layer (2,783 events with junction reads) [h4_summary.tsv, h4rob_summary.tsv, h4new_summary.tsv, logs/h4_adjustment_check.txt, logs/gse289499_depth_check.txt].

## Conclusion

For an RBP the model has never seen, a linear baseline matches the deep models in direction accuracy and beats them in rank correlation, and no model separates from its shuffled-RBP control by more than two points. The one class predicted well, the core spliceosome, is predicted mostly through an RBP-independent event prior. The measured responses themselves are informative: a splicing shift names its RBP across cell lines far above chance, and in tumours the one frequent, strong splicing driver tested (SF3B1 hotspot, uveal melanoma) retrieves its own class. This is the splicing-layer counterpart of the expression-layer result of Ahlmann-Eltze, Huber and Anders (Nature Methods 2025) that deep perturbation models do not consistently beat linear baselines on unseen perturbations.

## Reproduce

```
sbatch slurm/42_cnn_seq.sh      # 900-nt event sequences
sbatch slurm/44_ridge.sh        # ridge baseline, 20 (cell, fold) runs
sbatch slurm/43_train_cnn.sh    # v4 CNN, 20 runs
sbatch slurm/47_programme_prep.sh && sbatch slurm/48_tcga_coexpr.sh && sbatch --export=ALL,TAG=_final slurm/49_programme_model.sh   # programme-level model (in this order)
envs/pb_gpu_venv/bin/python 13_model/benchmark/03_benchmark_metrics.py
```
Outputs: `13_model/benchmark/master_table.tsv`, `master_by_rbp.tsv`, `paired_differences.tsv`, `difficulty_table.tsv`, `fig_benchmark.png`.
