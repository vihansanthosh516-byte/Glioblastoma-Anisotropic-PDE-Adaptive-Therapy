# Analysis plan v1 (internal preregistration)

Date written: 2026-10-04. Status: DRAFT until committed. After commit, change only by a dated amendment at the bottom.
Source: `C:\Users\vihan\Downloads\GBM_ISEF_PHD_Level_Masterplan.md` (§72, §73, §75, §107, §119).
Evidence rule: every number cites `output/` or `vault/PAPER_FINDINGS_LEDGER.md`.

## What was already seen before this plan (so it cannot be called blind)
- MU-Glioma core forecast (script 81): Dice 0.259 vs 0.233, n=133 (`output/forecast_labels_core/results.json`). The core target was picked after a volumes-only probe (ledger B2). This is development data. It is exploratory, not confirmatory.
- LUMIERE **volumes** (script 91): growth model lost to no-change (`output/lumiere_volume_forecast.json`). Already seen. Any LUMIERE volume claim is exploratory.
- LUMIERE **masks** (30 GB zip): NOT downloaded, NOT seen. This is the only clean confirmatory external test.
- rho > 0.02 cutoff: chosen after seeing data (ledger 2j). Post hoc.

## Research question
Can a patient-specific mechanistic model, updated from longitudinal MRI, forecast the next scan better than simple baselines, say when it can be trusted, and support treatment policies that stay good when the model is wrong? (§3)

Q1 Forecasting. Q2 External generalisation. Q3 Reliability. Q4 Decision-making (simulation only).

## Datasets and roles
| Role | Dataset | Notes |
|---|---|---|
| Development | MU-Glioma-Post | patient-held-out folds |
| External test | LUMIERE masks | frozen model only, no tuning |
| Stress test (optional) | UCSF longitudinal set or RHUH-GBM (n=40, `data/external/rhuh/`) | call it "distribution shift", not GBM confirmation |
| Molecular | multiomic-gbm, TCGA, CGGA | patient-held-out; not matched to MRI patients |
| Treatment | virtual cohorts + real growth rates | simulation; kill rates assumed |

## Primary endpoint (one only)
Per-patient paired difference in next-scan spatial overlap, delta_i = Dice(model_i) − Dice(persistence_i), on the contrast-enhancing core, MU development cohort, rolling-origin forecasts, patient as the unit.
- Estimate: mean and median delta, patient-level bootstrap 95% CI (10,000 resamples), % improved, % worse.
- Test: paired permutation test; Wilcoxon signed-rank as robustness.
- Success for H1: CI lower bound > 0 AND median delta > 0. Mean alone is not enough (the current median is +0.005, ledger B2).
- Confirmatory external endpoint: same delta on LUMIERE masks, frozen model. Direction and CI reported separately from development.

## Secondary endpoints (one family, Holm correction)
1. Volume error (symmetric percentage error)
2. HD95 / average surface distance
3. Centroid displacement
4. Interval coverage at 50/80/90/95% (calibration)
5. Growth-rate forecast error
6. Catastrophic-error rate (threshold fixed before running, below)
7. Selector vs fixed model at matched coverage

## Exploratory endpoints (no confirmatory language; BH-FDR if screened)
Survival links, molecular→PDE prior, treatment policy comparisons, RL vs MPC, DTI ablation, compartment forecasting, neural residual, value of information.

## Hypotheses (§119)
| ID | Test | Success condition |
|---|---|---|
| H1 | PDE vs persistence, MU, held-out patients | primary endpoint above |
| H2 | LUMIERE masks, frozen | same sign as H1, CI reported |
| H3 | anisotropic vs isotropic | external gain, CI > 0 |
| H4 | DTI vs non-DTI | external gain, CI > 0 (current: −0.001, p 0.97) |
| H5 | sequential update vs fixed fit | lower next-scan error |
| H6 | uncertainty vs error | monotone relation (Spearman, patient-level) |
| H7 | selective prediction | better at matched coverage |
| H8 | molecular prior | external forecast gain |
| H9 | robust MPC vs baselines | higher utility across uncertainty |
| H10 | RL vs robust MPC | RL beats it |
| H11 | probing | value of information > 0 |
| H12 | simulator shift | policy stays effective on generator B |

Not all are expected to be positive. Every result is logged, including negatives.

## Models compared (the ladder, §7–8)
Persistence; geometric expansion; last-rate / linear / exponential / Gompertz extrapolation; isotropic PDE; anisotropic PDE; DTI PDE; treatment-aware PDE; molecular-prior PDE; hybrid residual (last, only if external gain).

## Splits and leakage rules
- Patient is the independent unit everywhere (§76). Cells, voxels, and scans are not units.
- Folds assigned by patient ID before any fitting. Same split manifest for all models.
- Preprocessing, thresholds, feature selection, graph construction fit on training patients only.
- Forecast at scan t uses only scans ≤ t. Subgroups use only pre-forecast information.
- Tests that must pass: disjoint patient IDs, no duplicate scans, no temporal leakage (`tests/test_no_leakage.py`).

## Pre-set choices (fixed now; do not change after seeing external data)
- Horizon bins: short < 60 d, medium 60–120 d, long > 120 d. (Check these against the MU interval distribution before the first run; if the data make them empty, amend dated, before any model result.)
- Catastrophic error: Dice < 0.1 on a patient whose persistence Dice ≥ 0.1. Report the fraction for each model.
- Bootstrap: patients resampled, 10,000 draws, seed 20261004.
- Multiple testing: Holm for the secondary family; BH-FDR at 0.10 for exploratory screens.
- Subgroups: growth / stable / shrink defined from the previous two scans' volume change (cutoffs ±10%); early post-RT = within 12 weeks of radiotherapy end (RANO 2.0).
- Freeze: git tag `v1.0-development-frozen` before the LUMIERE mask run. Changes after the run are a new exploratory analysis, filed separately.

## Stopping and selection rules
- Model choice (selector, thresholds, hyperparameters) made on development folds only.
- External cohort is run once with the frozen config. Re-runs only for crashes, logged.
- If numerical verification fails (Phase 2), no biological result is reported until fixed.

## Wording rules
Say "simulation", "assumed", "post hoc", "exploratory", "retrospective". Do not say "clinical tool", "validated digital twin", "treatment recommendation". (§109–110, §154)

## Amendments
(none yet)
