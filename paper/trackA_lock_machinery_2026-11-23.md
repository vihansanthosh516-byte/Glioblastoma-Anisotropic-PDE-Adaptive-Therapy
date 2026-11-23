# Track A lock machinery (handoff section 5)

2026-11-23. Code `src/trackA_lock.py`, tests `tests/test_trackA_lock.py`. No real data used. Scratch-copy run: 24 of 24 tests passed (11 new lock tests + 13 earlier development tests), Python 3.14.

## What it does

| Handoff step | Implemented as |
|---|---|
| Score one independent unit | `regional_unit_scores` (1 / 0.5 tie / 0 per patient); `tmz_unit_errors` (absolute origin-mean error); `paired_differences` (larger favours candidate for both endpoints). Every enrolled unit must have a prediction, so fallback cases stay in the denominator. |
| Comparison family | `decide` takes all primary comparators; each must clear the gate. |
| Pass rule | One-sided paired t lower bound per comparator at alpha / k (Bonferroni, simultaneous); pass only if all exceed `min_gain`. Exact sign test with Holm adjustment is reported separately as a win/loss description, not the pass rule. |
| Size the rule | `power` simulates the actual `decide` rule with correlated comparators, ties and attrition; `units_for_power` finds the unit count for a target. |
| Runner | `preflight` (one role per origin, known roles) → `save_predictions` (same units across arms, no outcome keys, finite) → `sealed_evaluate` (verifies signed protocol, config and input hashes, fixed gates; requires saved predictions; writes an outcome-access receipt before reading outcomes; refuses reruns). All writes atomic and non-overwriting. |

Tests cover: ties, missing units, TMZ sign direction, Holm, every-comparator rule, power rising with n and staying low under no effect, unsigned protocol, changed input hash, unfilled gate, outcome access before predictions, mismatched arms, outcome-like keys, NaN, interrupted write, duplicate origin across roles, and rerun refusal.

## Illustrative sizing (NOT a protocol value)

Under an assumed paired-difference SD of 0.08, 3 comparators with correlation 0.5, 10% attrition, alpha 0.05, min gain 0 and 80% power: about 99 independent origins for a mean gain of 0.03 and about 35 for 0.05 (400 simulation draws each; Monte Carlo error a few units). These inputs are invented for the illustration. Real values must come from permitted development estimates.

## Not done

- No signed protocol exists for either endpoint. Gates, comparator list and power inputs are still open.
- Source-specific adapters for Park, Nam, Petti or GBM-Space data wait on the metadata replies.
- Environment and data-version capture in the receipt is limited to hashes and commit; a full environment lock file is still needed.
