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

## Amendment 1 (2026-10-04, before any model result is run under this plan)
Triggers: self-review (checked against `run_improved_aniso.py`, `src/81_label_correct_forecast.py`, `output/forecast_labels_core/`) and an outside review. Items below SUPERSEDE earlier text where they conflict. Existing development-cohort results (script 81) predate this plan and stay exploratory.

### A1.1 Unit of analysis and aggregation
- Patient is the only independent unit. A patient with several rolling-origin forecasts gets ONE value: delta_i = mean over that patient's eligible forecasts of (Dice_model - Dice_persistence), equal weight per forecast.
- All CIs: resample patients (10,000 draws), carrying all of a patient's forecasts together. Scans and voxels are never resampled.
- Per-forecast analyses (error vs interval, subgroups) use the cluster bootstrap by patient.
- Script 81 used only scan index 0 -> 1 (one pair per patient, 152 pairs, 133 scored). Phase 1 first counts patients with >= 3 scans. If fewer than 30 patients have >= 3 scans, rolling-origin is reported as secondary and the primary uses the first pair per patient.

### A1.2 Eligibility (primary population)
Source: `output/mu_glioma_cohort.json`. Include a patient if ALL hold:
1. >= 2 timepoints with non-null `day_from_diagnosis`; 2. interval > 0 days; 3. the target mask (labels 1 and 3, the core) is non-empty at the scan used as input; 4. brain mask and segmentation files load.
Exclusions are listed with a reason in the CONSORT table. Empty target at the scan being forecast is kept and scores Dice 0 for every arm (as in script 81). No patient is removed after seeing a result.

### A1.3 Failure rules (fixed now)
| Event | Rule |
|---|---|
| missing scan day / missing file | excluded, listed |
| interval < 14 d or > 365 d | kept in the secondary all-intervals analysis; excluded from the primary population; both reported |
| empty mask at 2 mm after resampling | excluded if input scan, Dice 0 if target scan |
| registration or tensor-atlas failure | patient excluded from every arm, listed |
| PDE solver NaN / non-finite / step-limit hit | that forecast scores Dice 0 for the PDE arm only, counted and reported; no re-fit by hand |
| fit lands on a parameter-grid edge | flagged per fold and reported; no manual change |

Nothing is repaired using outcome information.

### A1.4 Baselines (every arm uses the same crop, domain, 2 mm grid, 0.5 visibility, scorer)
Let V_t, V_{t-1} be volumes at the input scans and dt the gap in days. S(t) is the input mask.
- Persistence: forecast mask = S(t).
- Geometric expansion: forecast = S(t) dilated (distance transform inside the brain mask) until its volume equals V_hat. V_hat comes from information before the forecast only. Script 81's `uniform_dilation` uses the true scan-2 volume, which is an oracle. Keep it as an upper bound labelled "oracle", never as a baseline.
- Last-rate: g = ln(V_t / V_{t-1}) / (days between); V_hat = V_t * exp(g * dt). Needs >= 2 prior scans; for a first forecast it falls back to persistence, flagged.
- Linear: V_hat = V_t + (V_t - V_{t-1}) * dt / (days between), floored at 0.
- Exponential: same as last-rate (reported once).
- Gompertz: fit V(t) = K * exp(ln(V0/K) * exp(-a t)) on prior scans only, K fixed at the training-fold 95th percentile volume; needs >= 3 prior scans.
- PDE arms: below.

### A1.5 PDE specification (frozen as implemented in `run_improved_aniso.py` on 2026-10-04)
- Equation: du/dt = div(D grad u) + rho u (1 - u) - k(t) u, with k(t) = alpha * C_TMZ(t) + beta * R(t), alpha = 0.08 /day, beta = 0.03 /day per Gy/day (assumed, from `treatment_aware_pde.py`).
- Tensor arms: aniso (DTI atlas, sharpening r in {1, 10}), iso_same (atlas mean diffusivity), iso_homog (constant 1). D is scaled by d.
- Domain: brain mask AND atlas tissue, OR the scan-1 core. Boundary: zero flux on faces touching voxels outside the domain.
- Initial condition: u0 = 1 inside the input core mask, 0 elsewhere. Crop: bounding box of the core + 20 voxels.
- Discretization: finite volume, 2 mm voxels, explicit Euler, dt = min(0.5 d, 0.9 h^2 / (6 lambda_max), 0.05 / (rho_max + k_max)), u clamped to [0, 1].
- Visibility: forecast mask = u > 0.5.
- Fitting: no optimizer. Exhaustive grid rho in {0, 0.01, 0.03, 0.1} /day, d in {0.01, 0.03, 0.1, 0.3} mm^2/day, r in {1, 10}. Choice by 5-fold patient-level CV, best mean Dice on the training patients, scored out-of-fold. No regularization, no stopping rule.
- KNOWN PROBLEM: 4 of 5 folds chose rho = 0.1 (top of grid) and d = 0.01 (bottom of grid). Before the freeze, run a development-only grid-extension check (rho up to 0.3, d down to 0.003). If the best fit moves, the extended grid becomes the frozen grid and the primary result is re-run on development data. This is a development step, not external tuning.
- Solver verification (Phase 2) must pass before any biological number is reported.

### A1.6 H1 wording
"Broadly positive" requires BOTH mean delta CI lower bound > 0 AND median delta > 0 (patient level). One without the other is "mixed". This stops a few large wins driving the result. Current development numbers (mean +0.026, median +0.005) are not assumed to hold after the re-run.

### A1.7 H2 external (LUMIERE masks)
- Report LUMIERE mean and median delta and the 95% patient-level CI whatever the result.
- Strong replication: point estimate > 0 and CI wholly above 0. Weak: point estimate > 0, CI includes 0. Not replicated: point estimate <= 0. The loss of the growth model on LUMIERE volumes (ledger 2l) is a prior warning, not a reason to change the rule.
- No tuning on LUMIERE: not one parameter, threshold, or preprocessing step.

### A1.8 H3 anisotropy and the frozen set
All MU-fitted choices (grid, selection rule, aniso/iso arms, thresholds, preprocessing, DTI atlas build) are frozen from MU before any LUMIERE run. If LUMIERE has no tensors, H3 and H4 are untestable externally and are reported as "untested", not as null. Development data already show iso_homog (+0.028 vs no-change) at least as good as aniso (+0.026), and aniso vs iso_same is -0.001 (p 0.97): the gain is not from direction.

### A1.9 H4 DTI becomes exploratory
Moved to the exploratory list. The development result (-0.001, p 0.97) is reported as seen. No formal success condition.

### A1.10 H8 molecular prior (rewritten)
H8: a population-derived molecular prior (program scores derived in TCGA / atlas / CGGA, patient-held-out) shifts the prior distribution of rho for MRI patients. There is no patient-level link between the molecular and MRI cohorts. The test is whether the population prior changes held-out MRI forecast error versus a flat prior. It cannot show that a patient's own molecular state improves their forecast. If no matched molecular-imaging data are found, the result is labelled exploratory.

### A1.11 H9 utility (fixed before any policy run)
U = w1 * TTP/T - w2 * mean(M(t))/M(0) - w3 * cumulative dose / dose_MTD - w4 * resistant fraction at T.
Horizon T = 180 days; no discounting; TTP censored at T; M = total tumour mass in the simulator. Progression = mass above 1.2 x baseline for 14 days (set now; amend only before policy runs).
Proposed weights w = (1, 1, 0.5, 0.5), set by us, not by data. The weight sweep and Pareto frontier are reported; the primary conclusion is for the proposed weights only. w3 is a "treatment-burden proxy", not toxicity.

### A1.12 H10 and information
H10 (neutral): compare standard schedule, MTD, a simple rule, robust MPC, and RL under matched information and constraints. "RL adds nothing" is a valid result. Oracle RL appears only as a labelled upper bound.
Information, identical for every non-oracle controller at each step: noisy tumour mass (and a spatial mask for spatial controllers), treatment history, elapsed time, and a posterior over (rho, D, k) updated from those observations. No controller sees the true kill rate, true resistant fraction, or true rho. All share the same action set, dose cap, observation interval, and noise.

### A1.13 Subgroups (mathematical)
For a forecast from scan t with a previous scan t-1: c = (V_t - V_{t-1}) / V_{t-1}. shrink: c < -0.10; stable: -0.10 <= c <= +0.10; growth: c > +0.10. Uses only scans up to t. First forecasts have no previous scan: they form a "no-history" stratum and are NOT assigned to growth / stable / shrink. The "grew / shrank-or-same" split in script 81 uses scan 2 (the target), so it is outcome-defined and descriptive only.

### A1.14 Early post-RT
"Early post-RT" = forecast input scan within 84 days after the last radiation day in the cohort table. It is an analysis-defined subgroup motivated by RANO 2.0, not a RANO category. Missing radiation dates go to an "unknown timing" stratum.

### A1.15 Catastrophic error
Dice < 0.1 on a forecast whose persistence Dice >= 0.1. It is an analysis-defined research threshold, not a clinical failure threshold. Sensitivity to 0.05 and 0.2 is reported.

### A1.16 Multiple testing
Primary endpoint alone, no correction. Secondary family: Holm. Exploratory screens: BH-FDR at 0.10. Methods are not changed after results.

### A1.17 Freeze checklist (before the LUMIERE mask run)
Git tag `v1.0-development-frozen`; the exact commit hash written in `configs/development.yaml`; frozen: model configs, hyperparameters, thresholds, preprocessing, feature selection, subgroup and bin cutoffs, selector; analysis code frozen except bug fixes (each logged with a diff and a reason). The run manifest records seed, package versions, and data hash.

### A1.18 External reruns
One normal LUMIERE run. A rerun is allowed only for a crash or environment error, and the crash is logged. A scientific result never triggers a rerun. Any change after the run is a new exploratory analysis in a separate folder.

### A1.19 Numerical verification gate
If solver verification fails (manufactured solution, convergence, positivity, conservation), stop all biological analysis until fixed, then re-run affected development results.

### A1.20 Horizon bins (verified on data)
Intervals in script 81 (n=152): 3 to 1109 d, median 72. Bins: short < 60 d (59), medium 60-120 d (66), long > 120 d (27). The primary population drops intervals < 14 d or > 365 d (A1.3). The long-bin CI will be wide. Bins stay as set.
