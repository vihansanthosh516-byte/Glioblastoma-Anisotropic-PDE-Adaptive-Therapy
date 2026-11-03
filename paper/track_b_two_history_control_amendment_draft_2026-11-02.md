# Track B amendment draft: two-history non-PDE control (PROPOSED, 2026-11-02)

**Status: PROPOSED. Not adopted.** This draft authorizes no job, no test, no patient-data read and no scoring. It changes no code or protocol file. Duke stays NO-GO (`paper/track_b_duke_go_no_go_packet_2026-11-02.md`). LUMIERE patient scoring stays paused.

**Superseded for future consideration by r2:** [`track_b_two_history_control_amendment_r2_draft_2026-11-02.md`](track_b_two_history_control_amendment_r2_draft_2026-11-02.md). The original §4 uses a cross-outer-fold fallback that can leak held-out outcomes; its volume fallback is undefined; and it mixes Duke label rules with the LUMIERE runner. Do not adopt the original text.

Sources: rule and constants `src/fk_two_history_shape.py`; spec and screen `paper/track_b_two_history_control_spec_2026-11-02.md`; screen output `output/track_b_two_history_screen_20261102/screen.json`.

## 1. The question the primary contrast answers

Two arms get the same patients, the same folds, the same support S, the same target labels and the **identical** volume draws:
- **Candidate:** the PDE-ranked shape.
- **Control:** a two-history signed-distance extrapolation.

They differ only in **how voxels are ranked** to fill each volume draw. The contrast answers one question:

> At identical history-only volume forecasts, does the PDE-derived ranking give higher next-visit core Dice than a two-history geometric extrapolation?

It does **not** isolate a causal PDE effect. Why:
- The PDE posterior and the boundary model bring their own fitted parts.
- The control's settings are selected by a different rule.
- The labels are automated.

So the claim is predictive: one system outperformed a specified non-PDE shape rule under shared constraints.

## 2. Volume forecast (shared by every shape arm)

- **Source:** the candidate's training-only volume trend in `src/fk_hybrid_forecast.py` (`requested_volume_draws_voxels`). It is fitted on outer-training patients only. For inner selection (§3), it is refitted on inner-training patients only.
- **Inputs:** history masks h1 and h2 and their intervals, plus the realized origin-to-target interval. No target volume, target mask or outcome of the patient being forecast is used.
- **Identical draws:** for each patient, the set of volume draws (in voxels) is computed once. The same draws go to the candidate, the control, the last-mask signed-distance arm and the matched-volume ablation. Each arm matches each draw exactly by top-k inside S.
- **Estimand:** the next observed visit, conditional on its realized interval. This is not a fixed-day forecast.
- **Failure:** if the volume forecast fails for a patient, no shape arm gets draws. Every arm then uses the same frozen fallback, and the patient stays in the denominator. This rule must be written in the adopting amendment **before** the control is wired.

## 3. Nested, patient-grouped settings selection

- **Grid (fixed before any data):** SHRINK_TAU_DAYS {0, 7, 28} × MAX_FACTOR {1, 4, 8} × MAX_SD_DELTA_MM {10, 30}. That is 18 cells, plus f = 0 (the last-mask rule) as a 19th cell. The grid is chosen for coverage. It is not the screen's winner. The screen's best cell (tau 28, cap 1, clip 10) gets no special status.
- **Outer loop:** the protocol's outer folds, grouped by patient.
- **Inner loop:** inside each outer-training set only, patient-grouped inner folds (proposed: 5, or leave-one-patient-out if fewer than 10 patients). For each inner split:
  1. Refit the volume trend on the inner-training patients.
  2. Forecast the inner-validation patients.
  3. Score every cell by pre-boundary patient-macro Dice.
- **Selection rule:** pick the cell with the highest mean inner-validation Dice. Ties within 0.005 go to the most damped cell, in this order: f = 0, then the lowest cap, then the highest tau, then the lowest clip.
- **Never used for selection:** any outer-test patient's target, volume or outcome. Outer-test patients are scored once, with the cell chosen for their fold.
- **Boundary perturbation:** this is not part of selection. Selection uses pre-boundary scores only.

## 4. Selection stability and the small-sample rule

Report for each outer fold:
- the selected cell
- inner mean Dice for every cell
- the margin between the best cell and f = 0
- the inner patient count

Also report how many outer folds picked the same cell.

Stability rule (proposed; it must be fixed before any data):
- **Unstable** if any outer-training set has fewer than 10 patients, or if fewer than 3 of 5 outer folds select the same cell.
- **If unstable**, the control falls back to **f = 0 (last-mask signed-distance)** in every fold. The primary contrast is then relabeled "candidate vs last-mask signed-distance at identical volume draws". The two-history claim is dropped. The fallback and its reason are reported.
- No settings change after outer-test scores are seen.

## 5. Scores kept apart

| Score type | What it is | Use |
|---|---|---|
| End-to-end pre-boundary Dice | Shared history-only volume draws, no boundary noise | **Proposed primary** contrast (candidate vs control) |
| End-to-end boundary-perturbed Dice | Same draws plus each arm's boundary perturbation (shared seeds) | Secondary. The boundary covariance is PDE-calibrated, so this contrast is not neutral. |
| Oracle-volume shape diagnostic | Each shape arm given the target's true volume | Diagnostic only. Never a forecast score. Reported in a separate table. |

Pre-boundary and boundary-perturbed scores are never pooled.

## 6. What the synthetic screen does and does not show (preserved)

- The `-0.114` (growth) and `-0.153` (shrinkage) copy-last differences compare copy-last with no volume scaling against last-mask signed-distance at the **true** target volume. They are an oracle shape diagnostic, not a forecast gain.
- The noisy draws there are **target-centered synthetic noise**, not realistic forecasts.
- The control's default was not stronger than last-mask signed-distance in growth (-0.007) or shrinkage (-0.036) pre-boundary.
- The screen's fixed synthetic boundary covariance does not show how the PDE-calibrated boundary model behaves on patients.
- Nothing from the screen is carried into this protocol.

## 7. Multiplicity

The adopting amendment must state the family before data is read. Proposed:
- **One primary contrast:** candidate vs control, end-to-end pre-boundary.
- **Max-T family:** the four existing comparators (`COMPARATORS` in `src/fk_lumiere_runner.py`) plus the control, adjusted together.
- **Descriptive only:** everything else, including the boundary-perturbed results, the oracle diagnostics and the stability tables.

The amendment must also record whether earlier MU or LUMIERE ablation results were seen before this control was designed.

## 8. The limited claim the study could support

If the primary contrast passes:

> "On [cohort] (exploratory; n patients; automated source labels), at identical history-only volume forecasts, the PDE-ranked forecast had higher next-visit core Dice than a two-history geometric extrapolation whose settings were chosen by nested patient-grouped validation."

It could not support:
- a causal PDE effect
- agreement with true tumor
- clinical value
- superiority over MRI-intensity image AI
- generalization beyond the cohort
- a novelty claim

A pass is not promised. A null result is plausible.

## 9. Current blockers

- **Duke access:** pending, not granted. No written grant is on file (`paper/track_b_duke_source_ledger_2026-11-02.md`).
- **School approval:** user-reported, not passed. No written copy with date, approver and conditions is saved.
- **Duke eligible triplet count:** unknown (blocked). There is no authorized metadata.
- **Duke rights:** CC BY-NC-ND 4.0 allows making adaptations for noncommercial use but does not authorize sharing them. Whether forecasts, masks or weights count as adaptations, and how the custom Duke license treats them, is unknown.
- **Compute:** the projection is a scenario, not a bound (`paper/track_b_compute_audit_2026-11-02.md`). Nested selection adds control fits that are not counted yet. These are cheap compared with the PDE calls, but they are not measured.

## 10. Conditions before any patient run

1. Written Duke access grant and the applicable terms, saved in the repo.
2. Written school approval with date, approver and conditions, saved in the repo.
3. Authorized metadata-only eligibility count that meets a predeclared patient minimum.
4. A dated, cohort-specific amendment that adopts:
   - this section's rules, or revised ones
   - the fallback rules
   - the stability rule
   - the multiplicity family
   - the label map (core = source labels 1|3, cavity = 4)
   - the support and the estimand
5. The control wired into the runner, with synthetic tests that check identical draws and that no outer-test data reaches selection. These tests must pass in a scratch copy.
6. An updated compute projection that includes nested selection.
7. A pre-run seal (`src/fk_track_b_run_guard.py`) written at a clean commit, with its SHA-256 recorded in the amendment.
8. User authorization for the patient run.

Until all eight are met, no patient run is authorized.
