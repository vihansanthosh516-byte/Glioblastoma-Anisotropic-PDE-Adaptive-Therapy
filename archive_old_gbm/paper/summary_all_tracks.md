> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Summary of all tracks

Written 2026-10-18. Built only from `README_STATUS.md` (2026-10-18) and `paper/overview_for_review.md` (older, about 2026-10-05). No script was re-run. Numbers marked "overview" have no file cited there and were not re-checked against `output/` or the ledger. Where the two files disagree, the newer one (`README_STATUS.md`) wins.

Scope of this archived GBM summary: simulation/retrospective public-data analysis, not a complete activity audit. Nothing was tested in a lab or on live patients. Paper status was author-reported at that snapshot; current dated state is ../../START_HERE.md.

## Rules of the study
- Prospective protocol intends declaration before result access and one frozen evaluation; this does not certify every historical access. H-2/H-3 were later reused as development under Amendment 12 and are SPENT, never future locked tests. Negatives are reported (`analysis_plan_v1.md`, `vault/failed_hypotheses.md`).
- Freeze files with hashes: `configs/h2_frozen.yaml`, `h3_frozen.yaml`, `h3v_frozen.yaml`. Git tags: `h2-frozen`, `h3-frozen`, `h3v-frozen`.

## Track A: single cell to clinical (`src/01-41`)
Data: multiomic-gbm atlas (140,355 cells; 15,000 used; 21 patients; zones Core, Periphery, Healthy) and TCGA-GBM (n=518). Source for data sizes: overview.

| Result | Number | Source | Weakest point |
|---|---|---|---|
| C-GAT zone classifier, random cell split | 78.7% vs scVI 73.0% | overview | Leak. Same patient in train and test |
| Patient-ID-only model | 82.8% | overview | Beats the C-GAT, so the split was leaking |
| Leak-free classifier (gat_pca_clean, patient-level folds) | 0.638 | `output/cgat_leak_free.json` | Lower than the leaking number |
| Patient-fold baselines | cells-per-patient 64.6%; PCA32+LR 66.5% | `output/trackA_audit.json:C_baselines_patient_folds`, ledger 2o | Donor/composition confounding, not biology proof |
| Age and survival, TCGA | HR 1.031/year, p 8.88e-16 | overview | Age is already well known |
| Age and survival, CGGA pooled complete-case adjusted n=343, events=294 | HR 1.011, p 0.018 | `output/cgga_adjusted_age_mgmt.json:Y5_age_adjusted.terms.age` | Different covariate/complete-case model from pooled unadjusted and IDH-wt subgroup; observational association |
| Age, CGGA unadjusted check (n=374) | HR 1.009, Holm p 0.082 | overview | Pooled unadjusted n=374; not the IDH-wt primary subgroup (HR1.018, Holm p=.016). Adjusted analysis is a separate model, not an automatic replacement. |
| MGMT test | detects only HR 1.74 or more | overview | Underpowered |
| Energy landscape | Historical two-state pseudo-trajectory model | `output/spib_saddle_point_metrics.json`, ledger A2 | Internal numerical check, not biological-state proof |
| Virtual drug screen | in silico only | overview | No lab test |
| Link from Track A to Track B | Default-valued score 1.0; no measured link established | ledger 2ap/2m | No ID linkage; likely synthetic zone inputs lack exact provenance manifest |

## Track B: tumour growth and adaptive therapy (`src/42-48` and later)
Historical counts: MU 154 fitted/133 scored and UCSF DTI 62 require exact cohort/inclusion manifests; source count is not a universal denominator. Forecast primary pairs and endpoint-specific exclusions are in ledger 2p/2r/2ak.

| Result | Number | Source | Weakest point |
|---|---|---|---|
| Forecast on core target vs no-change | Dice +0.026, Holm p 0.0015 | `output/forecast_labels_core/results.json` | Target chosen post hoc. Loses on shrinking cores (-0.039, n=50) |
| Same forecast, overview version | 0.259 vs 0.233, CI 0.008-0.044, n=133, median gain +0.005 | overview | Median gain is tiny; target selected after volume probe, low absolute Dice, not directly comparable to later split/window |
| Registered old-solver PDE vs persistence (H1), GBM-only n=105, aniso r=1 | +0.0130, CI -0.0038 to +0.0292 | `output/pde_manifest/results.json`, ledger 2y | CI includes 0 |
| Same, all eligible n=133 | +0.0088 | same | CI includes 0 |
| Monotone sensitivity DTI direction (aniso r=1 vs iso_same, n=105) | -0.0004 | `output/pde_manifest_monotone/results.json`, ledger 2ak | Negative sensitivity contrast, not the registered primary result |
| Stratified forecast | AUC 0.46 | `output/forecast_stratified/` | Below chance |
| Growth-only forecast vs no-change, LUMIERE | worse (p 0.0003 HD-GLIO-AUTO; p 0.004 DeepBraTumIA) | `output/lumiere_volume_forecast.json` | Negative |
| H-3v volume forecast | dev null; former LUMIERE test worse; later spent | ledger 2ba | Negative |
| Inverse estimation | D at a bound in 147/154; rho at lower bound in 124/154 | `output/inverse_est_metrics.json` | D not identifiable |
| Adaptive vs MTD, 61 simulated patients | MTD ends with more resistant cells in 54/61 (7 tied, 0 reversed), E_MAX ratio 1000 | `output/emax_sweep_summary.json` | Assumed kill rates. VR 2.39 (larger final tumour). 10/61 progress earlier |
| Adaptive dose used | 32.36% of MTD | `output/emax_sweep_summary.json`, ledger 2n | Assumed treatment model |
| Progression earlier, growth rate > 0.02/day | 10/10 earlier; 0/51 below | overview | Post hoc cutoff; 7/10 at OFF=.50 in later controller sensitivity (2bn) |
| Kill-scale sweep (ratios 100-8000) | adaptive never selects more resistance (0/61 at all 7); final-mass cost 1.00 to 3.17; same progression pattern only at 500-2000 | `output/emax_sweep_summary.json` | Exploratory |
| Controller-threshold sweep (8 cells) | 0/61 more resistant; under 50% of MTD dose in every cell; pattern changes at OFF 0.50 and at ON 1.2 with OFF 0.80/0.90 | `vault/audit_2026-10-18/P4_threshold_sweep_results.csv`, ledger 2bn | Exploratory |
| Elongation from anisotropy | +0.178, p 4e-19 (Wilcoxon, n=103) | `output/fractal_aniso_vs_iso.json` | Simulation vs simulation. Secondary endpoint. Primary (fractal dimension) went the other way (-0.018) |
| Elongation, overview version | 1.21 vs 1.03 | overview | Simulation vs simulation |

## Track C: digital twin and RL (`src/49-65` and later)
Data: real fitted growth rates, assumed kill rates, synthetic patient sets.

| Result | Number | Source | Weakest point |
|---|---|---|---|
| RL vs Stupp, simulated cohort64 (n=20) | PPO wins 25% | `output/rl_kill_conditioned/evaluate.json`, ledger C1 | Simulated |
| RL on 21 real fitted-growth patients (real_test) | PPO wins 20/21 on day-90 volume; never on time to progression | same | Kill rates assumed. A one-line rule matches it |
| Paced dosing rule | Shared-assumed-kill fitted-input check: about +2.1 days (ledger 2j/C3); old randomized-kill ranges unverified here | `output/limit_checks.json`; exact cohort/arm keys required for `output/probe_paced_policies/evaluate.json` | Different kill assumptions/cohorts are not pooled; no measured patient gain |
| Probe-then-commit | 100% in stated noise-free model; 62% is cohort64 sigma=.10 rule arm, 3-seed mean (45/75/65%) | ledger C2 | Cohort/arm-specific simulation; other cohorts differ |
| Dual-drug rescue | came from the second drug, not adaptive control | overview | Not a control result |
| "Real patient" set | real growth rates, same assumed kill rates for all | overview | Not truly personal |
| Measured segmentation noise | floor 0.88 (middle tertile) vs 1.02 assumed; by tertile 1.03 / 0.88 / 0.35 | `output/seg_noise_floor.json` | Two automatic tools, not test-retest |
| Per-patient treatment effect | MGMT TTP p 0.60; resistance gain 0 of 6 | `output/mgmt_resistance/`, `output/plastic_resistance/results.json` | Negative |

Do not call this a "digital twin". There is no validated forecast and no two-way data flow.

## Historical recurrence study (PREDICT-GBM, H-2 to H-4)

Both H-2/H-3 former tests are SPENT under Amendment 12. Primary enhancing endpoints failed; positive secondary all-recurrence estimates are not fresh locked validation or clinical benefit. H-4 remains post-hoc development without an eligible new test.
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
| Released U-Net rebuilt | author-reported reproduction on eight dev patients; exact metric/hash artifact unverified here | script 125 |

Main limit: in 32 of 49 UCSF patients the 15 mm plan already covers at least 0.999 of enhancing tissue.
Retrain caveat: attempt 2 is a second launch of a "one-shot" run. Ledger 2bg discloses it.

## Historical dataset search for H-4

Rejection/access claims below are desk notes; dated search protocol and primary eligibility records are not linked. They are not confirmed current availability. See corrected current data-source notes before reuse.
- Rejected: CFB-GBM v2.0 (no pre-op MRI, one merged GTV). Burdenko-GBM-Progression (post-op baseline, restricted licence, no enhancing label).
- Desk-checked only, none confirmed: BraTS 2025, UCSD-PTGBM, Duke longitudinal (Zenodo 20548528).

## Open items
1. Find an H-4 locked test set (A13.4). If none, apply A13.7: H-4 stays "post hoc dev".
2. Write the paper. References are unchecked.
3. Check ChatGPT feedback against files; record in `paper/chatgpt_review_response.md`.
4. Historical 2026-10-17 dependency/storage inventory is superseded or unverified; current optional extras cover some named imports. Exact clean-install and storage measurements require checks; do not quote those old sizes as current.
5. Ledger 2bo is author-reported split runs totaling 101 without retained transcript. Current Git-backed scratch: 119 passed, 13 warnings, 27.74 s (`../../paper/codex_audit2/fix_progress.md`); test success is not scientific validation.
