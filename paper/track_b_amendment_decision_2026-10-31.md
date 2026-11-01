# Track B runtime amendment: decision memo (2026-10-31)

**Status: amendment A-2026-10-31-R ADOPTED by the user on 2026-11-01 (reuse only, full-fingerprint match, central-gap separate). The crop is NOT adopted. Adoption does not authorize scoring.** Wired: `add_learned_baselines_arm(unet_reuse_cache=...)`; the protocol wrapper shares one cache across the three prior runs and gives the central-gap refit its own. Earlier status line kept below for the record.

~~Status: draft for the user's decision. NOT ADOPTED.~~ No patient data was read or scored. The runner's active recipe is unchanged (`UNET_RECIPE`, `crop_voxels=None`). All timing below is synthetic and is not patient-result evidence.

Sources: `paper/track_b_runtime_pilot_2026-10-31.md`, `output/track_b_timing_pilot_20261031/timing_record.json`, `output/track_b_timing_pilot_20261031/crop_timing_record.json`.

This memo supersedes draft A-2026-10-31-U in the pilot report. That draft combined reuse and the 64³ crop. The replacement below, A-2026-10-31-R, is reuse only.

## 1. Recommendation: reuse the U-Net across the three prior-support runs, with a full fingerprint match

**Recommend: adopt reuse.** It saves U-Net CPU time and does not change any U-Net forecast when the fingerprint matches.

Verified identical across the primary, one-decade-narrower and one-decade-wider prior runs (code check of `run_mu_first_two_to_third_with_prior_sensitivity` in `src/fk_mu_patient_runner.py`; detail in the pilot report):

| Item | Same? | Evidence |
|---|---|---|
| Inputs (visit-1 core, visit-2 core, origin domain) | yes | same cohort records and domains |
| Targets (visit-3 core ∩ domain) | yes | same records |
| Folds | yes | same `cohort.fold_by_patient` |
| Preprocessing (padding) | yes | depends on masks only |
| Seeds | yes | 20261102 + fold + offset; fold counter restarts per runner call |
| Recipe | yes | `UNET_RECIPE`, not a runner option |
| Runtime environment | **checked by fingerprint, per fold** | see below |
| PDE prior bounds | not used by the U-Net | U-Net fit never reads the PDE nuisance |

**Fingerprint extended (code, this commit).** `src/fk_learned_baselines.py:unet_runtime_environment` records torch version, torch git build, CUDA build, numpy version, device (`cpu`), intra-op and inter-op thread counts, deterministic-algorithms flag during training (forced `True`) and the MKL-DNN flag. `unet_training_fingerprint` now hashes this environment together with patient IDs, input and target arrays, fold and recipe. Each fold also stores `runtime_environment` beside `training_fingerprint`.

Tests (scratch copy, Python 3.14, 50 passed across `tests/test_fk_learned_baselines.py`, `tests/test_fk_fixed_growth.py`, `tests/test_fk_mu_pipeline_components.py`):
- `test_unet_identical_across_prior_runs_but_not_gap_subset`: two prior runs give equal fingerprints and equal masks; a patient subset does not.
- `test_fingerprint_covers_runtime_environment`: a change in torch version, device, thread count or deterministic flag changes the fingerprint.

**Rule:** reuse a fold model only when the complete fingerprint matches. On any mismatch, retrain that fold. The central-gap refit has a different patient set, so it always trains its own models.

**Weakest point:** the fingerprint checks recorded settings. It does not prove bitwise-equal weights on other hardware (for example, a different CPU instruction set with the same torch build). Reuse inside one run on one machine avoids this. Reuse is now wired into the runner (2026-11-01).

## 2. The crop cannot be the U-Net novelty bar

**Keep the full-grid U-Net as the required comparator.** The 64³ crop is not adopted (user decision 2026-11-01).

Why:
- The crop limits where the U-Net can predict. Fixed growth, PDE-only, copy-last and the candidate keep the full support. So the crop can only lower the U-Net's Dice. That makes the U-Net an easier baseline to beat. It favours the candidate.
- Dropping patients with target outside the crop does not fix this. It changes the cohort, and the dropped patients are the ones with distant growth. It cannot replace the full-cohort primary comparison.
- A crop is fair only if every method uses the same predeclared region and the same scoring rule. That is a different estimand and needs its own amendment before any target mask is opened.

So a win over a cropped U-Net is not a win over the required U-Net baseline.

## 3. CPU-hour needs (synthetic; HISTORICAL MU 18-patient cohort)

**[Historical, MU 18-patient synthetic pilot on a 120×120×80 grid; superseded for LUMIERE by `paper/track_b_no_contact_dataset_plan_2026-10-31.md` §7a.]** Current LUMIERE projection: 677.6 CPU h with reuse, 1,124.1 h without; PDE 57.8 h per analysis is a lower bound (`output/track_b_timing_pilot_20261031/lumiere_timing_record.json`).

**Measured** (synthetic cohort, this PC, `crop_timing_record.json` unless noted):

| Item | Value |
|---|---|
| Full-grid U-Net step, 120×120×80 | 5.340 s mean (5.081–6.881), 40 steps |
| 64³ crop U-Net step | 1.735 s mean (1.545–2.920), 40 steps |
| Full-grid step, first pilot (`timing_record.json`) | 6.94 s mean (5.82–8.97), 12 steps |
| PDE full forecast, 1 patient (`timing_record.json`) | 3,128.6 s |
| PDE calibration-only forecast, 1 patient (`timing_record.json`) | 1,808.7 s |

Step time varied by about 25% between the two pilots.

**Projected** (not measured; 64,800 U-Net steps per analysis at the mean step; four analyses; PDE never reused):

| Scenario | Protocol-complete CPU hours |
|---|---|
| Frozen full grid, no reuse | 591.7 |
| Frozen full grid, with reuse (MU, historical) | 399.5 |
| 64³ crop, no reuse (not recommended as primary) | 332.2 |
| 64³ crop, with reuse (not recommended as primary) | 269.7 |

Every row includes the PDE at 207.3 h over four analyses (51.8 h × 4). **The PDE figure is a lower bound.** It counts 18 full and 72 calibration-only forecasts per analysis; inner cross-fitting may need more. With reuse, the central-gap U-Net is counted at full size.

Track B has no fixed time cap. These are resource needs, not a pass or a fail.

## 4. Count model zero-event correction

**Exact formula** (`src/fk_learned_baselines.py:fit_count_model`). With training counts y_i and gaps t_i (days) from outer-training patients only:

  mean_rate = (Σ y_i) / (Σ t_i)  if Σ y_i > 0
  mean_rate = 0.5 / (Σ t_i)      if Σ y_i = 0

The mean rate is the fallback rate and the reference rate. The Poisson IRLS start value uses log(max(Σ y, 0.5) / Σ t); this only sets the start point.

**Rationale.** With zero training events the rate would be 0. Any test patient with one or more new pieces would then get infinite Poisson deviance. Adding 0.5 to the zero total is a standard half-count continuity choice. It is training-only smoothing: it uses no test count and no test target.

**Status.** Set during the synthetic pilot, before any patient outcome was seen. It is frozen. Do not tune it against test counts. It is recorded in `paper/track_b_fixed_growth_baseline_spec_2026-10-31.md` §5.

**Weakest point:** 0.5 is a convention, not an estimate. If many folds have zero training events, count deviance will depend on it. Report the number of folds that used it.

## 5. Amendment A-2026-10-31-R (reuse only; ADOPTED 2026-11-01)

> **Track B amendment A-2026-10-31-R (draft).**
>
> (1) The fold-trained U-Net is trained once per outer fold in the primary-bounds analysis. The one-decade-narrower and one-decade-wider prior-support analyses reuse that fold's models unchanged, but only if the stored training fingerprint matches exactly. The fingerprint covers patient IDs, input and target arrays, fold, the full U-Net recipe, and the runtime environment (torch version and build, numpy version, device, thread counts, deterministic setting, MKL-DNN flag). On any mismatch the fold is retrained.
>
> (2) The central-gap refit trains its own U-Net models.
>
> (3) Nothing else changes. The U-Net stays full-grid (no crop), 300 epochs, three seeds, widths 16/32/64/128, Adam 1e-3, threshold 0.5. All required mask baselines (copy-last, fixed growth, PDE-only, fold-trained U-Net), all endpoints (patient-macro plain Dice on the full-grid future core, paired differences, count metrics) and the pass rule stay as frozen.
>
> Adopted on 2026-11-01, before any Track B target mask is opened.

**Effect.** No forecast and no endpoint changes. Only CPU time changes (MU historical 591.7 → 399.5 h; LUMIERE 1,124.1 → 677.6 h). Fairness is unchanged.

**Active:** reuse is wired with the fingerprint check and tested (synthetic only). Adoption date recorded.

## 6. Gates still closed

- Source-backed valid target support.
- MRI-only longitudinal registration QC.
- Cavity and label-zero meaning, from source documents (not from target voxels). Fail closed without it.
- Protocol freeze and hashes.
- School-use approval: user-reported; written copy not saved. (Duke access applies to the MU/Duke path only.)

The runner still refuses to score until the loader records the registration and target-support gates as passed. No patient mask has been read or scored.

## Decisions (2026-11-01)

1. A-2026-10-31-R adopted (reuse only). Crop not adopted. Full-grid U-Net required.
2. Patient voxel reads and scoring stay blocked until the protocol is frozen and the user approves.
