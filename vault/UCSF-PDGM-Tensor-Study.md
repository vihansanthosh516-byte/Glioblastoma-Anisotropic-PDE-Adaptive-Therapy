# UCSF-PDGM Three-Way Tensor Study (real patients, n = 62)

Written 2026-10-04. This study was missing from the vault. Numbers below were recomputed from the CSV files on disk.

**Question:** does a PDE driven by an anisotropic tensor predict the tumour mask better than the same PDE with an isotropic tensor?
**Metric:** Dice between predicted and real segmentation, paired per patient.
**Split:** 42 train / 20 test, seed 42. Source: `output/cv_split.json`.

| Mode | Dice aniso | Dice iso | Delta | Aniso wins (of 62) | Source |
|---|---|---|---|---|---|
| Standard DTI (patient's own eigenvalues) | 0.7843 (README) | 0.9464 (README) | -0.162 | **0** (train 0/42, test 0/20) | `output/results_standard.csv`, `cv_train_results.csv`, `cv_test_results.csv` |
| Atlas-enhanced DTI | **0.8229** | 0.7987 | +0.0242 | **12** | `output/results_atlas.csv` |
| Kurtosis-adjusted (DKI) | 0.1065 | 0.0275 | +0.0790 | **60** | `output/results_dki.csv` |

Train / test deltas (from `cv_*_results.csv`): Standard -0.162 / -0.162. Atlas +0.018 / +0.036. DKI +0.067 / +0.103.

## Corrections to the README
- README says Atlas aniso wins "~35 / 62". **The CSV says 12 / 62.** Use 12. (Train 8/42, test 4/20.)
- Standard DTI values 0.7843 / 0.9464 come from the README. The CSV for Standard (`results_standard.csv`) was not re-averaged in this audit. The train and test files give 0.7848 / 0.9468 and 0.7833 / 0.9457, which agree. Recompute before quoting 4 decimals.

## How to read it honestly
- Raw patient DTI makes anisotropic modelling **worse** than isotropic in every patient.
- The "wins" for Atlas and DKI come from heuristic eigenvalue scaling (lambda1 x1.5 etc.). They are tuned choices, not patient biology.
- DKI absolute Dice is 0.03-0.11. Both arms are poor. The win is in a regime where neither model matches the mask. Do not present it as a good prediction.
- The cohort is grade 2-3 IDH-mutant glioma, **not GBM**. State this.
- No external replication.

## Link to the forecast result
Script 81 (MU-Glioma-Post, GBM) finds the same story: DTI orientation adds nothing (aniso - iso_same = -0.001, p = 0.97). See [[Script-81-Label-Correct-Forecast]].

## Files
`run_improved_aniso.py` (driver, about 50 KB), `cross_validation.py`, `isef_figures.py`. Figures: `output/fig1_bar_chart.png`, `fig2_scatter.png`, `fig3_delta_boxplot.png` (not opened in this audit).
