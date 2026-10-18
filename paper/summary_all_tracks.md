> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Summary of all tracks

Written 2026-10-18. Built only from `README_STATUS.md` (2026-10-18) and `paper/overview_for_review.md` (older, about 2026-10-05). No script was re-run. Numbers marked "overview" have no file cited there and were not re-checked against `output/` or the ledger. Where the two files disagree, the newer one (`README_STATUS.md`) wins.

All work is simulation or retrospective analysis of public data. Nothing was tested in a lab or on live patients. The paper is not written yet.

## Rules of the study
- Each test is declared before data is read. A locked test set is used once. Negatives are reported (`analysis_plan_v1.md`, `vault/failed_hypotheses.md`).
- Freeze files with hashes: `configs/h2_frozen.yaml`, `h3_frozen.yaml`, `h3v_frozen.yaml`. Git tags: `h2-frozen`, `h3-frozen`, `h3v-frozen`.

## Track A: single cell to clinical (`src/01-41`)
Data: multiomic-gbm atlas (140,355 cells; 15,000 used; 21 patients; zones Core, Periphery, Healthy) and TCGA-GBM (n=518). Source for data sizes: overview.

| Result | Number | Source | Weakest point |
|---|---|---|---|
| C-GAT zone classifier, random cell split | 78.7% vs scVI 73.0% | overview | Leak. Same patient in train and test |
| Patient-ID-only model | 82.8% | overview | Beats the C-GAT, so the split was leaking |
| Leak-free classifier (gat_pca_clean, patient-level folds) | 0.638 | `output/cgat_leak_free.json` | Lower than the leaking number |
| Simple models, patient-level split | 60-66% | overview | Not re-checked |
| Age and survival, TCGA | HR 1.031/year, p 9e-16 | overview | Age is already well known |
| Age and survival, CGGA adjusted | HR 1.011, p 0.018 | `output/cgga_adjusted_age_mgmt.json` | About 1/3 of the TCGA effect size |
| Age, CGGA unadjusted check (n=374) | HR 1.009, Holm p 0.08 | overview | Older number. The adjusted HR above is newer |
| MGMT test | detects only HR 1.74 or more | overview | Underpowered |
| Energy landscape | 2 metastable states | overview | No test or file cited |
| Virtual drug screen | in silico only | overview | No lab test |
| Link from Track A to Track B | score was 1.0 for every patient | overview | IDs do not match. File is probably synthetic |

## Track B: tumour growth and adaptive therapy (`src/42-48` and later)
Data: MU-Glioma-Post MRI (154 fitted, 133 scored) and UCSF-PDGM diffusion MRI (62 patients).

| Result | Number | Source | Weakest point |
|---|---|---|---|
| Forecast on core target vs no-change | Dice +0.026, Holm p 0.0015 | `output/forecast_labels_core/results.json` | Target chosen post hoc. Loses on shrinking cores (-0.039, n=50) |
| Same forecast, overview version | 0.259 vs 0.233, CI 0.008-0.044, n=133, median gain +0.005 | overview | Median gain is tiny |
| PDE beats persistence (H1), GBM-only n=105, aniso r=1 | +0.0130, CI -0.0038 to +0.0292 | `output/pde_manifest/results.json`, ledger 2y | CI includes 0 |
| Same, all eligible n=133 | +0.0088 | same | CI includes 0 |
| DTI direction (aniso r=1 vs iso_same, n=105) | -0.0004 | `output/pde_manifest_monotone/results.json`, ledger 2ak | DTI adds nothing |
| Stratified forecast | AUC 0.46 | `output/forecast_stratified/` | Below chance |
| Growth-only forecast vs no-change, LUMIERE | worse (p 0.0003 HD-GLIO-AUTO; p 0.004 DeepBraTumIA) | `output/lumiere_volume_forecast.json` | Negative |
| H-3v volume forecast | dev null; locked LUMIERE worse than persistence | ledger 2ba | Negative |
| Inverse estimation | D at a bound in 147/154; rho at lower bound in 124/154 | `output/inverse_est_metrics.json` | D not identifiable |
| Adaptive vs MTD, 61 simulated patients | MTD ends with more resistant cells in 54/61 (7 tied, 0 reversed), E_MAX ratio 1000 | `output/emax_sweep_summary.json` | Assumed kill rates. VR 2.39 (larger final tumour). 10/61 progress earlier |
| Adaptive dose used | 32% of MTD | overview | Not re-checked |
| Progression earlier, growth rate > 0.02/day | 10/10 earlier; 0/51 below | overview | Cutoff chosen after seeing data |
| Kill-scale sweep (ratios 100-8000) | adaptive never selects more resistance (0/61 at all 7); final-mass cost 1.00 to 3.17; same progression pattern only at 500-2000 | `output/emax_sweep_summary.json` | Exploratory |
| Controller-threshold sweep (8 cells) | 0/61 more resistant; under 50% of MTD dose in every cell; pattern changes at OFF 0.50 and at ON 1.2 with OFF 0.80/0.90 | README_STATUS (file not named) | Exploratory |
| Elongation from anisotropy | +0.178, p 4e-19 (Wilcoxon, n=103) | `output/fractal_aniso_vs_iso.json` | Simulation vs simulation. Secondary endpoint. Primary (fractal dimension) went the other way (-0.018) |
| Elongation, overview version | 1.21 vs 1.03 | overview | Simulation vs simulation |

## Track C: digital twin and RL (`src/49-65` and later)
Data: real fitted growth rates, assumed kill rates, synthetic patient sets.

| Result | Number | Source | Weakest point |
|---|---|---|---|
| RL vs Stupp, simulated cohort64 (n=20) | PPO wins 25% | `output/rl_kill_conditioned/evaluate.json`, ledger C1 | Simulated |
| RL on 21 real fitted-growth patients (real_test) | PPO wins 20/21 on day-90 volume; never on time to progression | same | Kill rates assumed. A one-line rule matches it |
| "Paced" dosing rule, time to progression | +48 to +53 days (synthetic); +8 to +9 days (91 MU patients); +16 to +19 days (34-38 LUMIERE patients) | overview | Random kill rates. Not re-checked |
| Probe-then-commit | 100% no noise; 62% with realistic noise | overview | Not re-checked |
| Dual-drug rescue | came from the second drug, not adaptive control | overview | Not a control result |
| "Real patient" set | real growth rates, same assumed kill rates for all | overview | Not truly personal |
| Measured segmentation noise | floor 0.88 (middle tertile) vs 1.02 assumed; by tertile 1.03 / 0.88 / 0.35 | `output/seg_noise_floor.json` | Two automatic tools, not test-retest |
| Per-patient treatment effect | MGMT TTP p 0.60; resistance gain 0 of 6 | `output/mgmt_resistance/`, `output/plastic_resistance/results.json` | Negative |

Do not call this a "digital twin". There is no validated forecast and no two-way data flow.

## Main study: recurrence target (PREDICT-GBM, H-2 to H-4)
| Step | Result | Source |
|---|---|---|
| H-2, all-recurrence coverage | +0.081 (CI +0.065 to +0.097), n=101. Secondary. Positive | `output/predictgbm/final/test_results.json` |
| H-2, enhancing recurrence | +0.0088 (CI -0.0115 to +0.0307), n=98. Primary. Null | same |
| H-3 (UCSF), all-recurrence | +0.085 (CI +0.059 to +0.115), n=59. Secondary. Positive | `output/predictgbm/h3_final/ucsf_test_results.json`, ledger 2bd |
| H-3, enhancing recurrence | +0.0044 (CI -0.0182 to +0.0369), n=49. Primary. Null | ledger 2bd |
| H-3 at 90% volume | -0.0053. Non-inferiority at -0.01 not shown | ledger 2bd |
| H-3 dev arms (137 TUM patients) | C1 +0.030, C3 +0.027, C2 +0.030 (interim) | README_STATUS |
| H-4 step 1 (plan rules, pooled dev) | enh +0.0223, all +0.0869. Post hoc | ledger 2be |
| H-4 step 2 (retrain) | 636 samples, 15 epochs, seed 20261009. Loss 0.7994 to 0.7164. SHA-256 `b9846e1f...`. No test run | `kaggle_run/out_h4_retrain/`, ledger 2bg |
| H-4 step 3 (new locked test) | Blocked. No eligible dataset | README_STATUS 5a.1 |
| Released U-Net rebuilt | reproduces predictions on 8 dev patients | script 125 |

Main limit: in 32 of 49 UCSF patients the 15 mm plan already covers at least 0.999 of enhancing tissue.
Retrain caveat: attempt 2 is a second launch of a "one-shot" run. Ledger 2bg discloses it.

## Dataset search for H-4
- Rejected: CFB-GBM v2.0 (no pre-op MRI, one merged GTV). Burdenko-GBM-Progression (post-op baseline, restricted licence, no enhancing label).
- Desk-checked only, none confirmed: BraTS 2025, UCSD-PTGBM, Duke longitudinal (Zenodo 20548528).

## Open items
1. Find an H-4 locked test set (A13.4). If none, apply A13.7: H-4 stays "post hoc dev".
2. Write the paper. References are unchecked.
3. Check ChatGPT feedback against files; record in `paper/chatgpt_review_response.md`.
4. Repo cleanup: `pyproject.toml` lacks SALib, nibabel, dipy, scikit-learn, h5py. Untracked run kits: 2.0 GB untracked, 4.87 GB in `kaggle_run` (2026-10-17 audit).
5. Tests: 101 collected, 101 passed on a copy with partial data (ledger 2bo).
