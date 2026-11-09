# Related work for the PerturbBridge forward model (compiled 2026-11-09)

Every entry below was checked against a web search on 2026-11-09 for title, venue and link. Author lists marked
"(authors not confirmed)" were not visible in the search results and should be checked on the linked page before citing.

## A. Same task family: predict the splicing response to RBP knockdown from binding / sequence (ENCODE data)

1. Van Nostrand EL et al. A large-scale binding and functional map of human RNA-binding proteins. Nature 583, 711 (2020).
   https://www.nature.com/articles/s41586-020-2077-3 — the data we use (472 KD RNA-seq profiles / 263 RBPs; 223 eCLIP / 150 RBPs). Descriptive integration, no held-out-RBP prediction.
2. Van Nostrand EL et al. Principles of RNA processing from analysis of enhanced CLIP maps for 150 RNA binding proteins. Genome Biology (2020).
   https://www.encodeproject.org/publications/a93bb897-2c2c-49c1-944a-a4bbc8374c33/ — companion eCLIP analysis.
3. Raghav Y, Paul A, Anderson R, Karthyk S, Iturralde A, Vyas J, Dy J, Jones BC, Castaldi PJ, Platig J. Explainable machine learning reveals an RBP regulatory logic of exon skipping. bioRxiv (2026), under review.
   https://www.biorxiv.org/content/10.64898/2026.05.29.728731v1 ; code https://github.com/PlatigLab/ENCODE-RNA-Binding-Protein-Network-Modeling — eCLIP + shRNA KD RNA-seq in 2 cell lines, predicts PSI of skipped exons, XGBoost + SHAP; generalisation tested on unseen EVENTS, not unseen RBPs.
4. McGurk MP, McWatters DC, Burge CB. KATMAP infers splicing factor activity and regulatory targets from knockdown data. Nature Biotechnology (2025), doi 10.1038/s41587-025-02881-9; preprint bioRxiv 10.1101/2024.06.25.600605.
   https://www.biorxiv.org/content/10.1101/2024.06.25.600605v2.full — interpretable regression of splicing change on binding change, applied to HepG2/K562 KD data; infers activity per factor, not a transfer to unseen factors.
5. Gupta K, Yang C, McCue K, Bastani O, Sharp PA, Burge CB, Solar-Lezama A. Improved modeling of RNA-binding protein motifs in an interpretable neural model of RNA splicing. Genome Biology 25, 23 (2024).
   https://link.springer.com/article/10.1186/s13059-023-03162-x — 79 RBP in-vitro affinity models inside an interpretable splicing model.
6. Jha A, Gazzara MR, Barash Y. Integrative deep models for alternative splicing. Bioinformatics 33(14), i274–i282 (2017).
   https://pmc.ncbi.nlm.nih.gov/articles/PMC5870723 — RNA-seq + CLIP-seq, predicts tissue and KD/overexpression splicing changes for RBFOX/CELF/MBNL in mouse; the closest older precedent for "binding -> KD response".
7. Functional characterization of splicing regulatory elements. bioRxiv 10.1101/2021.05.14.444228 (2021) (authors not confirmed).
   https://www.biorxiv.org/content/10.1101/2021.05.14.444228.full.pdf — ENCODE eCLIP + KD, "sensitive vs robust" exons.
8. Large-scale evaluation of the ability of RNA-binding proteins to activate exon inclusion. Nature Biotechnology (2023/2024) (authors not confirmed).
   https://www.nature.com/articles/s41587-023-02014-0 — tethering screen, experimental ground truth for RBP splicing activity.
9. PARNET: a CLIP-seq-based foundation model for RNA sequence representation learning. bioRxiv (2026) (authors not confirmed).
   https://www.biorxiv.org/content/10.64898/2026.08.08.743506v1.full — 150 RBPs / 223 eCLIP; validated with KD RNA-seq "RNA maps"; binding model, not a KD-response predictor.

## B. Same evaluation problem, different layer: predict the response to an UNSEEN genetic perturbation (expression)

10. Roohani Y, Huang K, Leskovec J. Predicting transcriptional outcomes of novel multigene perturbations with GEARS. Nature Biotechnology (2023), doi 10.1038/s41587-023-01905-6.
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11180609/
11. Cui H et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. Nature Methods (2024), doi 10.1038/s41592-024-02201-0.
12. Hao M et al. Large-scale foundation model on single-cell transcriptomics (scFoundation). Nature Methods 21, 1481–1491 (2024), doi 10.1038/s41592-024-02305-7.
13. Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. Nature Methods 22, 1657–1661 (2025), doi 10.1038/s41592-025-02772-6.
    https://pmc.ncbi.nlm.nih.gov/articles/PMC12328236/ ; code https://github.com/const-ae/linear_perturbation_prediction-Paper
    Verified quote: "None of the deep learning models was able to consistently outperform the mean prediction or the linear model."
14. Kernfeld EM, Yang Y, Weinstock JS, Battle AJ, Cahan P. A systematic comparison of computational methods for expression forecasting. bioRxiv 10.1101/2023.07.28.551039 (2023, updated 2024-10-01; preprint only as of search).
    https://pubmed.ncbi.nlm.nih.gov/37577640/ ; benchmark site https://ekernf01.github.io/perturbation-benchmarks
    Verified quote (author blog): "Judging by most evaluation metrics, on most datasets, the mean and median perform better than most methods."
15. Csendes G, Sanz G, Szalay KZ, Szalai B. Benchmarking foundation cell models for post-perturbation RNA-seq prediction. BMC Genomics (2025), doi 10.1186/s12864-025-11600-2.
    Verified quote: "even the simplest baseline model—taking the mean of training examples—outperformed scGPT and scFoundation"; "current Perturb-Seq benchmark datasets exhibit low perturbation-specific variance".
16. Miller HE, Mejia GM, Leblanc FJA et al. Deep learning perturbation models can outperform baselines on calibrated metrics. Nature Biotechnology (2026), doi 10.1038/s41587-026-03307-w; preprint bioRxiv 10.1101/2025.10.20.683304.
    https://www.nature.com/articles/s41587-026-03307-w — the counter-argument: with calibrated metrics deep models do beat uninformative baselines.

## C. Sequence-to-splicing models (different task: no perturbation; useful as feature sources, not as score comparators)

17. Jaganathan K et al. Predicting splicing from primary sequence with deep learning (SpliceAI). Cell (2019), doi 10.1016/j.cell.2018.12.015. https://github.com/illumina/SpliceAI
18. Zeng T, Li YI. Predicting RNA splicing from DNA sequence using Pangolin. Genome Biology 23, 103 (2022). https://pmc.ncbi.nlm.nih.gov/articles/PMC9022248
19. Horlacher M, Wagner N, Moyon L, Kuret K, Goedert N, Salvatore M, Ule J, Gagneur J, Winther O, Marsico A. Towards in silico CLIP-seq: predicting protein-RNA interaction via sequence-to-signal learning (RBPNet). Genome Biology 24, 180 (2023), doi 10.1186/s13059-023-03015-7. https://pmc.ncbi.nlm.nih.gov/articles/PMC10403857

## What none of A–C does
No entry evaluates prediction of the knockdown splicing response for RBPs held out of training, with a shuffled-identity control and a
same-RBP independent-experiment ceiling. That is the gap PerturbBridge fills; its negative/limited result matches group B's conclusion.
