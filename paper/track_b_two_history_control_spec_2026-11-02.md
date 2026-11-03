# Track B two-history non-PDE shape control: development spec

**Status:** development specification for `src/fk_two_history_shape.py`, synthetic data only. Not a protocol amendment, not wired into `src/fk_lumiere_runner.py`, no patient data read. No learned residual exists or is claimed. The constants below are proposed development defaults; they are PROPOSED (not adopted) and become fixed only if a dated cohort-specific amendment adopts them before target scoring.

## Rule

With h1 (earlier history mask), h2 (origin mask), history gap `g = day(h2) - day(h1)`, forecast interval `T`, and signed-distance maps SD (positive inside, physical mm via the grid spacing):

`SD_f = SD_h2 + f * clip(SD_h2 - SD_h1, -30, +30)`, `f = min( g/(g+7) * T/g , 4 )`.

For each candidate volume draw k (the candidate's `requested_volume_draws_voxels`, same draws), the comparator selects the k highest-`SD_f` voxels inside the same support S (`volume_match_concentration`, stable flat-index ties), so each draw's volume matches exactly. Draws are combined as the candidate scorer does: mean over draws, inclusion probability >= 0.5, within S.

## Proposed (not adopted) constants and edge handling

| Item | Value / behaviour |
|---|---|
| `MAX_FACTOR` | 4.0 cap on f |
| `MIN_HISTORY_GAP_DAYS` | 1.0; gap below this (including 0) gives f = 0, i.e. the existing last-mask rule |
| `SHRINK_TAU_DAYS` | 7.0; weight g/(g+7) shrinks extrapolation when the history gap is short |
| `MAX_SD_DELTA_MM` | 30.0 clip on the SD difference before scaling |
| h1 empty or full | f = 0 (last-mask rule), flagged `h1_unusable_fallback` |
| h2 empty or full, h2 outside S, nonpositive T | error (same as the existing last-mask arm) |
| Shrinking history | SD difference is negative, so the shape contracts; volume still equals the draw |
| h1 == h2 | f*delta = 0, identical to the last-mask ranking (tested against `signed_distance_volume_baseline_draws`) |

The constants were chosen as plausible regularisers, not tuned; they were not fitted to any data. Their sensitivity has not been studied. A value could be revisited only on development cases before an amendment.

**Strength of the control is unverified.** A development sensitivity grid over tau, cap and clip, compared with the copy-last arm and the last-mask signed-distance arm (`signed_distance_volume_baseline_draws`), is required before adoption. Until it exists, a PDE win over this control may only mean the control is weak.

**f reduces to T/(g+7).** Since `g/(g+7) * T/g = T/(g+7)`, f is the forecast interval over the shrunk history gap, capped at 4. The history gap cancels in the numerator, so the control extrapolates by interval length, not by measured speed alone.

**Tie-break limit.** Equal SD_f values are broken by flat index (grid order). This is arbitrary and can place the volume boundary differently for equal-distance voxels.

**ForecastFailure fallback.** The fallback for the control (for example the last mask) must be defined in the amendment before wiring; none is defined here.

**Family and prior results.** The control is not in the max-T comparator family (`COMPARATORS` in the runner). The amendment must define the family it enters and must record whether earlier MU/LUMIERE ablation results were seen before the control and its constants were chosen.

## Boundary-perturbation matching: finding

**Exact matching of the displacement realisations is not possible.** In `correlated_boundary_perturbation` the displacement is `L(anchors(mask)) @ z(seed)`: the anchors are surface voxels of the volume-matched mask being perturbed, and the Cholesky factor `L` is built from those anchor positions. The perturbation then moves that mask's own signed distance and re-thresholds to the target count. The PDE-candidate's and the comparator's volume-matched masks differ, so they have different anchors, different `L`, and different displacements even with the same seed. (The perturbation acts on the rasterised matched-volume mask, not on the PDE field, but the mask-specific anchors are what break exactness.) Also, the existing no-PDE arm uses a separate seed stream (`boundary_seed.spawn(2)[0]`), not the candidate's per-draw seeds.

What is implemented:
1. **Best common-random-number version:** `candidate_boundary_seeds(seed, n)` reproduces the candidate's per-draw seeds (`SeedSequence(seed).spawn(2)[1].spawn(n)`; a test checks the derivation), so the standard-normal vector `z` is shared per draw. When the anchor count is equal (n_anchors capped by surface voxels), draws are correlated but not identical realisations. `n_anchors` must be passed the same as the candidate.
2. **Deterministic pre-boundary sensitivity:** `apply_boundary=False` returns the matched-volume masks before any perturbation. Compare against the candidate's pre-boundary quantity (`pre_boundary_probability`, sized masks before perturbation), both at the same volume draws. This isolates shape ranking from boundary noise. **Recommended: make the deterministic pre-boundary matched-volume contrast the primary or a co-primary contrast**, because the boundary covariance was calibrated on the PDE arm, so a boundary-perturbed contrast favours the PDE arm by construction. Report the boundary-perturbed contrast beside it.

Not tested here: the shared-seed path against a real candidate run (the PDE solver was not run); only the seed derivation and the comparator's internal behaviour are tested.

## What candidate-vs-control can and cannot show

- It can show whether, on a frozen cohort and endpoint, the PDE-derived spatial ranking gives different next-visit Dice (agreement with source-provided labels) than a history-aware, matched-volume, non-PDE ranking, under shared volume draws, support, interval and boundary model.
- It cannot show that the PDE is causal tumor biology, that the model is a validated digital twin, or clinical benefit. Both arms are predictive systems under shared constraints (same support, same volume draws), so a difference is attributable to shape ranking only to the extent the other pieces are matched; boundary realisations are only partly matched (above).
- A null or negative contrast is plausible: the comparator uses the same volume information, and the SD extrapolation is a simple, strong geometric prior. Report it unchanged.
- With only two history masks the comparator's velocity estimate is noisy; the cap and shrinkage are untuned.
- Patient, not voxel or visit, remains the inference unit; the cohort-specific primary-contrast promotion still needs a new dated amendment before target scoring.

## Tests

`tests/test_fk_two_history_shape.py` (12 tests): exact volume per draw, mask inside support, empty h1/h2, tiny and zero gap, factor cap and shrinkage, shrinking history, seed determinism and shared seeds, seed derivation, reduction to the last-mask rule when h1 == h2, difference from the last-mask arm when history moves, pre-boundary path, input validation, consensus mask. Run in a scratch copy only, never the working tree.
