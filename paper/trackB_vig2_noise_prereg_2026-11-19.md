# Track B: noise-aware level-reversion (VIG-3 candidates) — written before any run (2026-11-19)

Follows `paper/trackB_vig2_dev_result_2026-11-19.md`. Development only; no test is spent. Same six cohorts, same leave-one-cohort-out folds, same bootstrap (20,000, seed 20261119).

## Why this family, and its limit
Finding: level-reversion helped on cohorts with noisy input sizes (CFB tool labels, WAW, LUMIERE) and hurt on UCSF (accurate labels): H1 was +0.083 worse than persistence there. The family was chosen **after** seeing these results, so all six folds are partly informed. Only a future untouched cohort is clean.

I tried to use a cohort-level noise estimate (scatter of the earlier change over short gaps). It cannot be built: UCSF has no earlier scan for any pair, UCSD 14 of 55, WAW 14 of 48, so the estimate exists for three cohorts at most and could not be learned. Dropped before running anything.

## Candidates (scan-level, learnable, usable on a new cohort)
Base is H1 of the previous pre-registration: `y_in + g(dt) + (1 - exp(-a·dt))·(K - y_in) - h`, absolute loss. Changes:
- **K bounded** in every fold: `K = lo + (hi - lo)·sigmoid(k)` with lo = median and hi = maximum of `y_in` of the training cohorts of that fold. (H1's K went to 15.3, outside the data.)
- **N1 (size-dependent pull):** `a = a0 · sigmoid((y0 - y_in)/w)`. Small input cores (relatively noisier in log units) are pulled more. Parameters a0 > 0, y0, w > 0.
- **N2 (N1 plus history):** `a = a0 · sigmoid((y0 - y_in)/w) · exp(-b·has_history)`, b ≥ 0. A second reading (an earlier scan) should lower the pull.
- **N0 (control):** H1 with bounded K only, to show the effect of the bound alone.

References: persistence, VIG (`fit_vig`, unchanged), H1 (previous run, refitted by this script), Gompertz-L1.

## Selection rule (fixed now)
A candidate (N1 or N2) is selected only if ALL hold:
1. cohort-equal mean difference to persistence over six cohorts is lower than VIG's and lower than H1's;
2. its mean over the four earlier cohorts (MU, UCSD, UCSF, LUMIERE) is not worse than VIG's;
3. no cohort is worse than persistence by more than +0.03;
4. it beats persistence in at least 5 of 6 cohorts.
If none passes: nothing is frozen and the result is "no noise-aware gain". The lowest cohort-equal mean among those that pass wins. Nothing is added or retuned after the run. Starts and bounds are those in `src/track_b_vig2_noise_dev.py`.

## What a pass would mean
Only a candidate to freeze before a future untouched cohort. Not an improvement claim.
