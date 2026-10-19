> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# UCSF-PDGM historical three-way model comparison (62 retained result rows)

Written2026-10-04; source tables checked2026-10-19 without running the scientific driver. All six named split/result files exist in root output (audit absence claim incorrect). Means/win counts below are read-only summaries of retained CSVs, not independent replication.

**Question:** does a PDE driven by an anisotropic tensor predict the tumour mask better than the same PDE with an isotropic tensor?
**Metric:** Dice between predicted and real segmentation, paired per patient.
**Split:** 42 train / 20 test, seed 42. Source: `output/cv_split.json`.

| Mode | Dice aniso | Dice iso | Delta | Aniso wins (of 62) | Source |
|---|---|---|---|---|---|
| Standard DTI (patient's own eigenvalues) | 0.784316 (CSV) | 0.946408 (CSV) | -0.162 | **0** (train 0/42, test 0/20) | `output/results_standard.csv`, `cv_train_results.csv`, `cv_test_results.csv` |
| Atlas-enhanced DTI | **0.8229** | 0.7987 | +0.0242 | **12** | `output/results_atlas.csv` |
| Kurtosis-adjusted (DKI) | 0.1065 | 0.0275 | +0.0790 | **60** | `output/results_dki.csv` |

Train / test deltas (from `cv_*_results.csv`): Standard -0.162 / -0.162. Atlas +0.018 / +0.036. DKI +0.067 / +0.103.

## Corrections to the README
- README says Atlas aniso wins "~35 / 62". **The CSV says 12 / 62.** Use 12. (Train 8/42, test 4/20.)
- Standard result CSV mean aniso0.78431594/iso0.94640834, delta-0.16209240; all62 deltas are negative. Atlas12/62 wins (train8/42,test4/20); DKI60/62 (train40/42,test20/20). cv_split.json n_train42,n_test20,random_seed42 agrees. These checks do not establish diagnoses, original split timing or absence of model selection.

## How to read it honestly
- Raw patient DTI makes anisotropic modelling **worse** than isotropic in every patient.
- The "wins" for Atlas and DKI come from heuristic eigenvalue scaling (lambda1 x1.5 etc.). They are tuned choices, not patient biology.
- DKI absolute Dice is 0.03-0.11. Both arms are poor. The win is in a regime where neither model matches the mask. Do not present it as a good prediction.
- The old grade2–3/IDH-mutant cohort assertion needs an exact eligibility/clinical manifest; result CSVs alone do not establish it. Do not generalize this as GBM validation.
- No external replication is documented by this result bundle.

## Link to the forecast result
Script 81 (MU-Glioma-Post, GBM) is a different retrospective internal OOF study and does not demonstrate added DTI gain (aniso - iso_same = -0.001, p = 0.97). See [[Script-81-Label-Correct-Forecast]].

## Files
`run_improved_aniso.py` (driver, about 50 KB), `cross_validation.py`, `isef_figures.py`. Figures are historical references, not visually verified here; exact filename/byte provenance requires checking.

Exact source columns: results_{standard,atlas,dki}.csv dsc_aniso,dsc_iso,delta (62 rows each); cv_{train,test}_results.csv mode/delta (126/60 rows across three modes). Heuristic scaling is assumed, not a measured atlas or kurtosis correction. Underlying driver files exist at root, contrary to the audit search result; no patient-data driver or cross-validation job was rerun. No external replication or clinical treatment validation is established.
