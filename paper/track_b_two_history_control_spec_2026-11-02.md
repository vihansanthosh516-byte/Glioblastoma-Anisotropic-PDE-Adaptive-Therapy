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

## 2026-11-02 synthetic sensitivity screen (development only)

Script `src/fk_two_history_sensitivity_screen.py`, test `tests/test_fk_two_history_sensitivity_screen.py`. Output: `output/track_b_two_history_screen_20261102/screen.json` (unedited copy of the scratch run, 81 cases, 175 s, `py -3.14`). Defaults in `fk_two_history_shape.py` were not changed.

**Design.** 48^3 grid, 2 mm, fixed ellipsoidal support S. Three synthetic truth regimes (growth with anisotropic drift, shrinkage, plateau) x history 4/8/16 weeks x forecast 4/8/16 weeks x 3 seeds. Truth h1, h2, target come from a known synthetic ellipsoid process plus a persistent irregularity field; the target is used only to score and to supply the TRUE volume draw. Every arm gets the same 4 volume draws (true volume plus 3 log-normal noisy draws, sd 0.15). Arms: `copy_last` (h2 unchanged, not volume matched; no boundary field), `last_mask_sd` (f=0, volume matched; identical to a volume-matched copy-last, so there is no separate arm), `control_default` (tau 7, cap 4, clip 30) and an 18-cell grid (tau 0/7/28 x cap 1/4/8 x clip 10/30). Fields: `dice_*_pre` (apply_boundary=False quantity) and `dice_*_boundary` (shared candidate seeds, fixed synthetic covariance sigma 2 mm, tau 0.5 mm, ell 8 mm) are kept separate. Paired difference is arm minus `last_mask_sd`, per case.

**Results (mean Dice over cases; paired diff vs `last_mask_sd` in brackets; screen.json `by_regime`, true-volume draw, pre-boundary / boundary-perturbed).**

| Regime (27 cases) | last_mask_sd pre | copy_last pre | control_default pre | control_default boundary | best grid cell (tau 28, cap 1, clip 10) pre / boundary |
|---|---|---|---|---|---|
| growth | 0.892 | 0.778 (-0.114) | 0.885 (-0.007) | 0.856 (-0.011) | 0.903 (+0.011) / 0.874 (+0.007) |
| shrinkage | 0.797 | 0.643 (-0.153) | 0.761 (-0.036) | 0.711 (-0.051) | 0.803 (+0.006) / 0.763 (+0.001) |
| plateau | 0.997 | 0.999 (+0.001) | 0.997 (0.000) | 0.942 (+0.006) | 0.997 (0.000) / 0.943 (+0.007) |

(`last_mask_sd` means are `control_default` mean minus its paired diff; noisy-draw fields in screen.json show the same pattern.)

**Where the control helps or hurts.**
- Growth: `control_default` is about equal to `last_mask_sd` on average (-0.007 pre). It helps at long forecast with long or medium history (hist 16w / fc 16w +0.038; hist 8w / fc 16w +0.021) and hurts with short history and long forecast (hist 4w / fc 16w -0.043; hist 4w / fc 8w -0.059). Per-cell n is 3 cases, so these are noisy.
- Shrinkage: the default is worse than the last-mask rule in most cells (mean -0.036 pre, -0.051 boundary). Worst cell hist 16w / fc 16w -0.236, where extrapolating a contraction overshoots because the synthetic tumour hits its radius floor. Small gains at fc 4w with hist 16w pre-boundary (+0.033) disappear with the boundary (-0.022).
- Plateau: no difference pre-boundary (target nearly equals h2, Dice about 0.997 for all arms; this regime is nearly trivial and says little). Boundary-perturbed differences are small and mixed.
- Settings: the strongest shrinkage of extrapolation (tau 28, cap 1, clip 10) was best in every regime, and weak shrinkage with a high cap (tau 0, cap 8) was worst. The default is not the best cell, so in this synthetic family the control only matches or beats the last-mask rule when it extrapolates little.
- `copy_last` is clearly worse than volume-matched arms in growth and shrinkage, so most of the control's value over `copy_last` is volume matching, not extrapolation.

**Limits.** Synthetic ellipsoids with smooth irregularity are far simpler than real GBM change; the truth process, the shrinkage radius floor, noise level and boundary covariance are my choices and drive the results. 3 seeds per cell; per-cell SEs are in screen.json and are large. The boundary covariance is fixed and synthetic, not the PDE-calibrated one. Plateau is nearly trivial. No resection, edema, treatment or MRI-segmentation noise. The screen shows the control is not a uniformly strong baseline (it can lose to the last-mask rule), nothing more.

**Statement.** Synthetic results do not prove this control is a strong patient baseline. The settings are not tuned from this screen for any patient run. If a dated amendment adopted values taken from this screen (for example tau 28, cap 1, clip 10), those values would be post hoc relative to this screen, and that would have to be recorded; adopting the untouched defaults is not made safer by the screen, since they were shown not to be best. The primary/co-primary status of this control stays PROPOSED until a dated amendment defines it.

### Interpretation correction (2026-11-02, after `0c08d70`)

- **Oracle shape diagnostic, not a forecast gain.** The `-0.114` (growth) and `-0.153` (shrinkage) copy-last differences compare **unscaled** copy-last with `last_mask_sd` given the **target's true future volume**. They measure shape at an oracle volume. They are not end-to-end forecast gains, and no arm here forecasts volume.
- **Noisy draws are target-centered.** The 3 log-normal draws (sd 0.15) are centered on the target volume. Call them **target-centered synthetic noise**. They are not realistic predicted volumes: a real history-only forecast has bias as well as variance.
- **Boundary model.** The boundary-perturbed fields use a fixed synthetic covariance (sigma 2 mm, tau 0.5 mm, ell 8 mm). They do not show how the PDE-calibrated boundary model will behave on patients. Pre-boundary and boundary-perturbed comparisons stay separate and must not be pooled.
- **Status.** The control stays PROPOSED. Its default was not stronger than `last_mask_sd` in growth (-0.007) or shrinkage (-0.036) pre-boundary. The best grid cell (tau 28, cap 1, clip 10) was selected on this same screen. It is **not adopted**, and it must not be carried into a patient test.

**Before any future test, choose one path (decision for the user):**
1. **Descriptive diagnostic (recommended).** Report the control next to the other arms. It is not a primary or co-primary contrast. This needs no settings choice and adds no post-hoc risk.
2. **Development-fold selection.** A dated amendment fixes the settings grid and the selection rule first. Then settings are selected only on development training folds, never on test or target patients. Only after that may the control be a primary or co-primary contrast.

**If another synthetic screen runs:** predeclare a history-only volume forecasting process with realistic bias and variance (not target-centered). Give every shape arm the identical volume draws. Report every setting. Do not select a winner for the patient test from it.

No patient run is authorized. Duke remains NO-GO.
