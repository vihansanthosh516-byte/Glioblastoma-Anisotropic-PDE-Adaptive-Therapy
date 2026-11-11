# Track C and integration follow-up (2026-11-10)

No policy utility was computed and no policies were compared. Confirmatory evaluation stays stopped.

## Review findings checked

| Item | Result | Evidence |
|---|---|---|
| C-g GE-betting bisection monotonicity | **Supported empirically, not proven.** Grid scan (4,001 means) on 200 synthetic samples (n = 7, 14, 30, 100, 207; beta, near-constant, Bernoulli, narrow-uniform data): 0 accepted means outside the bisection interval, 0 accept sets with holes. Locked as `tests/test_trackc_ge_betting_shape.py` (24 samples, 801-point grid). Bets depend on the wealth path, so a proof is still needed in the independent statistical review | `src/trackc_ge_betting.py:64-171` |
| Cross-arm common random numbers in 211 | **Verified.** Every evaluation arm is scored with the same `eval_seed_key` (`211_...:1082`); noise is drawn before `policy.reset` inside the simulator (review C-j) | `211_...:794-805, 1082` |
| C-f signal-model normality | **Verified, narrower than stated.** Generator and likelihood use the same `f0_to_resistance_score` mapping (`trackc_integrated_signal_model.py:130-147`; generator line 57 of the function). The conditional likelihood given f0 is exact. Only the claim that the signals are marginally standard normal fails when a world's f0 law differs from the reference envelope; this affects the (q, rho) identifiability argument, not filter weighting | same file |

## Code added

- `src/trackc_integrated_driver.py`: reads protocol settings from the config at run time (closes the "config not loaded" part of gap C8 for the driver path); checks the V3 prior; builds collision-free world schedules from separate scenario and policy seed blocks; verifies every world's arm plan has the primary pair (finite-set minimax vs matched d = 0) and the controls (molecular-only, MRI-only, no-signal); checks completeness and order of outputs.
  - `mode="evaluation"` refuses while the config lists open gates or if the precision target is not 14 days. On the current config it always refuses.
  - `mode="mechanics"` refuses the real policy benchmark (also through `functools.partial` or `wraps`), so only stubs run.
- Single-signal control arms (`finite_set_minimax_molecular_only_mpc`, `finite_set_minimax_mri_only_mpc`) via a new `channels` option on `SignalAmbiguityBeliefMPC`; enabled only with `include_single_signal_controls=True`, which the world runner sets. Default behavior unchanged.
- Tests: `tests/test_trackc_integrated_driver.py` (9 tests). Scratch run of the Track C and driver files: 127 passed, 1 skipped, 3 expected failures (C9-C11).

## What still blocks evaluation (not relaxed)

1. Interval method: GE-betting is a candidate; independent statistical review is not done.
2. Precision: the 14-day half-width is not reachable with the reserved block (gap C4). The config audit gives 27,907 worlds per family for the Hoeffding bound; that is a count for that bound, not a universal minimum. The earlier 207-world figure only clears one worst-case lower-bound criterion and is not sufficient for the target.
3. Seed table and registry reservation for evaluation worlds; runtime cap; clean-install check.
4. U7 design choices (81 vs 135 models; 112-day forecast length).

These cannot be cleared by code in this repository. Items 1 and 2 need an outside reviewer and either the user's authorization of a large runtime or an explicit amendment of the target.

**Decision 2026-11-10 (user):** evaluation stays stopped; target and success rule unchanged.
