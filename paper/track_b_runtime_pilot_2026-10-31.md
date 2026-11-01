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
