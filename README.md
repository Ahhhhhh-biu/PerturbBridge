# PerturbBridge v2

Predicting splicing responses to unseen RNA-binding-protein (RBP) perturbations, and inverse-mapping
splicing shifts in premalignant / at-risk tissue onto RBP perturbations.

Design notes and status log are kept locally (Chinese); the English summary of methods and results is in `15_manuscript/MANUSCRIPT_DRAFT.md` and `RESULTS.md`.

| Directory | Content |
|---|---|
| `00_access_requests/` | Controlled-access applications (dbGaP, HTAN, EGA): research-use statements, checklists, status |
| `01_encode_metadata/` | ENCODE experiment/file tables (`experiments.tsv`, `files.tsv`), summary, download manifest |
| `02_encode_download/` | ENCODE RSEM transcript/gene quantifications (GENCODE V29) and eCLIP peak BEDs |
| `03_suppa_psi/` | SUPPA2 event definitions (GENCODE v29) and per-experiment PSI |
| `04_h2_check/` | H2 go/no-go: perturbation-specific variance of ΔPSI vs Δlog2FC |
| `05_precancer_cohorts/` | GSE274628 (BRCA-carrier fallopian tube), GSE166044 (pre-diagnosis breast), GSE69429 (STIC exon array): run tables, sample tables, STAR + rMATS outputs, cell-composition scores |
| `06_features/` | eCLIP event features (coordinate-keyed), RBP annotation (ESM2 embeddings, GenePT, motif 5-mer vectors, Pfam), resources (phastCons, MaxEntScan) |
| `07_novelty_notes/` | Novelty check and verified related-work list (`related_work_2026-11-09.md`) |
| `08_junction_psi/` | Junction-level (rMATS) cross-technology test on 12 shRNA/CRISPR pairs that decided the readout layer |
| `09_bam_layer/` | ENCODE BAM manifest, knockdown/control pairs, rMATS prep/post outputs for all 876 pairs |
| `10_event_universe/` | Coordinate-keyed event universe (546,917 events), event features, eCLIP features, per-run and training tables, 5-mer counts |
| `11_splits/` | Run tiers, paralog-grouped RBP folds (BLASTP), leakage report |
| `12_baselines/` | Measurement ceilings, shared-control artefact, control-centring, simple baselines |
| `13_model/` | Forward model (`01_train.py`), aggregation, calibration, ceilings, clean target, stratified report, expression layer (`10`–`12`), results tables and figures |
| `14_inverse/` | RBP response dictionary, H3 retrieval (ENCODE, cross-technology control, TCGA via Snaptron + cBioPortal), H4 precancer cohorts, HeyA8 and TMA validation, figure |
| `15_manuscript/` | Manuscript draft (`MANUSCRIPT_DRAFT.md`) |
| `envs/` | `pb_venv` (pip: pandas/scipy/statsmodels/SUPPA deps), `SUPPA/` (git clone), SRA-toolkit settings |
| `slurm/` | SLURM scripts (absolute paths) |
| `logs/` | SLURM logs; `JOBS.txt` = job id ledger |

Environment:
```bash
PB=/well/blagden/users/eic789/PerturbBridge
$PB/envs/pb_venv/bin/python          # project Python
$PB/envs/pb_venv/bin/python $PB/envs/SUPPA/suppa.py   # SUPPA2
module load SRA-Toolkit/3.0.10-gompi-2023a STAR/2.7.11b-GCC-12.3.0 rMATS-turbo/4.3.0-gfbf-2023a BEDTools/2.31.0-GCC-12.3.0
```

Status log: `STATUS.md`.
