# Track B runtime pilot (synthetic only, 2026-10-31)

**Synthetic timing only. No patient data was read. This is not patient-result evidence.** No epoch, grid, seed or model setting was changed because of this pilot.

Record: `output/track_b_timing_pilot_20261031/timing_record.json` (created 2026-11-01T01:56 UTC). Code: `src/fk_track_b_timing_pilot.py`. Run in the scratch copy with Python 3.14.4, torch 2.14.0+cpu, 8 CPU threads, Windows 11, Intel64 Family 6 Model 189.

## Setup

- Synthetic cohort: 18 patients (the frozen MU development cohort size). Grid 120 × 120 × 78 at 2 mm isotropic (SRI24 240 × 240 × 155 at 1 mm, resampled). Ellipsoid domain. Each patient has three nested ball-shaped tumor cores, with a 76–104-day target gap.
- Outer-training sizes assumed 14, 14, 14, 15, 15 (folds of 4, 4, 4, 3, 3).
- Each arm runs in its own process, so peak memory is per arm (Windows peak working set).
- U-Net: 12 timed training steps after 2 warm-up steps, with the frozen recipe (widths 16/32/64/128, Adam 1e-3, flips and axial rotations, Dice + BCE, deterministic algorithms). Padded input 120 × 120 × 80; 1,402,833 parameters. **A complete fold at 300 epochs × 3 seeds was not run**: the warm-up predicted about 26 hours for one fold.
- PDE hybrid: one full forecast and one calibration-only forecast for one synthetic patient, with the frozen settings (n_coarse 21, n_refine 11, n_draws 200, n_anchors 256, log-uniform prior on the primary bounds).

## Results by arm

| Arm | Measured | Peak memory |
|---|---|---|
| U-Net training step | mean 6.94 s, min 5.82 s, max 8.97 s over 12 steps (warm-up steps 21.5 s, 10.0 s) | 2.18 GB |
| U-Net prediction, 3 seeds, 1 patient | 5.17 s | (same process) |
| PDE hybrid, full forecast, 1 patient | 3,128.6 s (52.1 min) | 1.17 GB |
| PDE hybrid, calibration-only, 1 patient | 1,808.7 s (30.1 min) | (same process) |
| Fixed growth: radius selection on 15 patients | 15.0 s | 0.19 GB |
| Fixed growth: one forecast | 0.07 s | (same process) |
| Count model: new pieces for 18 patients + fit | 1.53 s + 0.001 s | 0.36 GB |

An earlier separate warm-up gave U-Net steps of 4.6–6.7 s, so the step time varies by about ±25%.

## Projection (cautious)

Assumptions: step time is constant across epochs, seeds and folds; nothing else runs on the machine; U-Net cost = (sum of outer-training sizes = 72) × 300 epochs × 3 seeds = 64,800 steps; PDE cost = 18 full forecasts + 72 calibration-only forecasts. Inner cross-fitting may call more PDE forecasts than this, so the PDE figure is a **lower bound**.

| Item | Hours |
|---|---|
| U-Net, one outer fold (15 training patients) | 26.0 |
| U-Net, one analysis (5 folds), mean step | 125.0 |
| U-Net, one analysis, slowest step | 161.4 |
| PDE hybrid, one analysis (lower bound) | 51.8 |
| **Total, one analysis** | **176.8 (≈ 7.4 days)** |
| **Protocol-complete run (four nested analyses: three prior supports + central-gap refit)** | **707.1 (≈ 29.5 days)** |

## Decision

The full protocol-complete run is **not practical** on this PC (about 29 days of continuous CPU). Following the stop rule, no patient-data run should start under the current settings. The settings are unchanged. Patient data may not go to Kaggle (only Track C may use Kaggle), so a cloud GPU is not an option without new permission.

## Proposed dated amendment (2026-10-31, NOT adopted)

This proposal needs the user's approval before it applies. Adopt it, if at all, before any target mask is opened.

1. **Reuse U-Net fold models across the three prior-support analyses.** The U-Net does not use the ρ–D prior, so it gives the same models in all three prior runs (same patients, folds, seeds and recipe). Train it once per outer fold for the primary cohort, and once more for the central-gap refit. This cuts U-Net training from four runs to two. It does not change any U-Net result.
2. **History-only U-Net region of interest.** Crop each patient's three input channels and target to a fixed 64 × 64 × 64-voxel (128 mm) box centred on the visit-2 core's centroid. Only history sets the box. Target core outside the box counts as a miss, which is reported. A U-Net step should then cost roughly 64³ / (120 · 120 · 80) ≈ 23% of now (unmeasured; voxel scaling assumed). Keep 300 epochs, 3 seeds, widths and optimizer.
3. **Rerun this synthetic pilot after 1 and 2**, and freeze the run only if the protocol-complete projection is below a limit the user sets.

PDE cost (about 52 hours per analysis, lower bound) is a separate problem. Options for the user: run the three prior-sensitivity analyses after the primary result is frozen; or keep the current settings and accept about 9 days for PDE across four analyses. This pilot does not choose between them.

## Limits

- Synthetic balls are simpler than real masks. PDE time may depend on tumor size and posterior shape.
- One U-Net process with 8 threads; the record does not capture other machine load.
- The PDE count of forecasts per analysis is an assumption; the true number may be higher.

---

## Second pilot: U-Net reuse check and measured 64³ crop (2026-10-31)

**Still synthetic only. No patient data was read. The crop and the reuse are PROPOSED, NOT ADOPTED. The frozen recipe is unchanged** (`UNET_RECIPE.crop_voxels = None`).

Record: `output/track_b_timing_pilot_20261031/crop_timing_record.json` (created 2026-11-01T02:48 UTC). Same PC, same scratch copy and same torch build as above. Command: `python -m src.fk_track_b_timing_pilot --crop-pilot`.

### U-Net reuse across the three prior-support runs: inputs verified identical

`run_mu_first_two_to_third_with_prior_sensitivity` (`src/fk_mu_patient_runner.py`) calls the runner three times with the same `cohort` and the same `runner_options`. Only `rho_bounds` and `D_bounds` change between the calls.

| Item | Same across the three prior runs? | Why |
|---|---|---|
| Training inputs (visit-1 core, visit-2 core, origin domain) | yes | same cohort records and `origin_domain_by_patient` |
| Training targets (visit-3 core ∩ domain) | yes | same records |
| Fold membership | yes | same `cohort.fold_by_patient`; folds run in order 0–4 |
| Preprocessing (padding, crop if any) | yes | depends only on the masks |
| Seeds | yes | each runner call builds a fresh wrapper, so the fold counter restarts at 0; seed = 20261102 + fold + offset |
| Recipe | yes | `UNET_RECIPE`, not a runner option |
| PDE prior bounds | **not used by the U-Net** | the U-Net fit never reads the base nuisance |

The central-gap refit uses a smaller patient set, so its inputs are different. It stays separate.

**Code check.** Each fold now records `unet_training_fingerprint`. This is a SHA-256 over the patient IDs, the input and target arrays, the fold and the recipe (`src/fk_learned_baselines.py:unet_training_fingerprint`). The test `test_unet_identical_across_prior_runs_but_not_gap_subset` runs the arm twice with different base nuisances. The fingerprints and the U-Net masks match exactly. A patient subset gives a different fingerprint.

**Limit.** On real hardware, identical models also need the same torch build and thread count. The fingerprint does not prove that.

Reuse is **not wired in**. If reuse is adopted, reuse a model only when its fingerprint matches.

### Measured timing (40 timed steps each, same session)

| Arm | Padded input | Step mean (min–max) | Predict, 3 seeds, 1 patient | Peak memory |
|---|---|---|---|---|
| Full grid, rerun | 120 × 120 × 80 | 5.340 s (5.081–6.881) | 2.88 s | 2.20 GB |
| 64³ history crop | 64 × 64 × 64 | 1.735 s (1.545–2.920) | 1.24 s | 0.95 GB |

**Measured crop speed-up: 5.340 / 1.735 = 3.08× per step.** This is lower than the 4.4× that voxel count alone suggested.

The full-grid step in this rerun (5.34 s) was faster than in the first pilot (6.94 s). Step time varied by about 25% between runs.

### Projection (not measured)

Assumptions:
- Each analysis has 72 training pairs × 300 epochs × 3 seeds = 64,800 U-Net steps, at the mean step time.
- PDE hours are the first pilot's lower bound (51.8 h per analysis), times four. The PDE is never reused.
- With reuse, the central-gap U-Net is counted at full size, so that figure is an upper figure.
- Prediction time is not added.

| Scenario | U-Net h per analysis | Protocol-complete CPU hours (4 analyses) |
|---|---|---|
| Frozen full grid, no reuse | 96.1 (max step 123.9) | 591.7 |
| Frozen full grid, with reuse | 96.1 | 399.5 |
| 64³ crop, no reuse | 31.2 (max step 52.6) | 332.2 |
| 64³ crop, with reuse | 31.2 | 269.7 |

With crop and reuse, the U-Net share is 62.4 h. **The PDE lower bound (207.3 h over four analyses) is now the main cost.**

Track B has no fixed time cap. These figures are resource needs, not a pass or a fail.

### Draft amendment wording (NOT ADOPTED; do not apply until the user approves)

> **Track B amendment A-2026-10-31-U (draft).**
>
> (1) The fold-trained U-Net is trained once per outer fold for the primary cohort. It is reused unchanged in the one-decade-narrower and one-decade-wider prior-support analyses, but only if its recorded training fingerprint matches. The central-gap refit trains its own models.
>
> (2) U-Net inputs and training targets are cropped to a 64 × 64 × 64-voxel (128 mm) cube. The cube is centred on the rounded centroid of the visit-2 core and shifted to lie inside the grid. Only the visit-2 core sets the cube. Predictions are mapped back to the full grid, with zeros outside the cube. Target core outside the cube stays in the target and counts as a miss; its voxel count is reported per patient.
>
> All other U-Net settings stay as frozen: widths 16/32/64/128, Adam 1e-3, 300 epochs, three seeds, 0.5 threshold.
>
> Adopted on [date], before any Track B target mask is opened.

**How the crop changes the estimand and fairness.**
- **Estimand:** it does not change. It stays plain patient-macro Dice of the full-grid tumor core, with all arms scored on the same full target.
- **U-Net handicap:** the U-Net loses the ability to predict tumor more than about 64 mm from the visit-2 centroid. Distant new tumor counts as a miss for the U-Net only. Fixed growth, PDE and persistence keep the full grid.
- **Who it favours:** the crop can only lower the U-Net's Dice, which makes the U-Net an easier baseline to beat. So it favours the candidate.
- **What to report:** per-patient target voxels outside the crop. If any patient has them, report the U-Net comparison with and without those patients as a sensitivity check. Do not claim a win over the U-Net that depends on those patients.
- **Reuse:** it does not change any U-Net forecast, so it does not affect fairness.

## Remaining gates (unchanged)

These are all still open:
- Source-backed valid support.
- MRI-only registration QC.
- Cavity meaning, from source documents.
- Duke access and school approval.

No patient mask has been read or scored. The PDE runtime (≥ 51.8 h per analysis) needs a separate decision.
