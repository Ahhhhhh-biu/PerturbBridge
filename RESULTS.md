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
envs/pb_gpu_venv/bin/python 13_model/benchmark/03_benchmark_metrics.py
```
Outputs: `13_model/benchmark/master_table.tsv`, `master_by_rbp.tsv`, `paired_differences.tsv`, `difficulty_table.tsv`, `fig_benchmark.png`.
