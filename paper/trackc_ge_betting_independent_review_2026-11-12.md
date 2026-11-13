# Independent review: GE-betting code vs the paper (2026-11-12)

Reviewer: an independent fresh-context AI reviewer (not a human statistician). It was told to be strict, did not write the code, and ran checks only in a temporary folder. This file records its report and the project's answer to each note. The reviewer returned the report as a message; it is copied here in short form.

Code: `src/trackc_ge_betting.py`. Paper: Martinez-Taboada and Ramdas (2026), arXiv:2608.21694v2.

**Correction.** `paper/trackC_five_decisions_evidence_2026-11-12.md` (D5) said the project's GE-betting is the Waudby-Smith and Ramdas (JRSSB 2024) "hedged capital" interval. That is wrong. The code follows Martinez-Taboada and Ramdas (2026). ConfSeq implements the 2024 hedged-capital method, so the ConfSeq comparison is a cross-check against a second valid method, not a check for equal numbers (see `output/trackc_interval_crosscheck/result.json`).

## Element check (reviewer)

| Element | Verdict |
|---|---|
| Variance path start 0.25 and prior mean 0.5; update rule (paper eq 1.3) | MATCH (numeric difference 0.0) |
| psi(p) = phi(Phi^-1(p)) / p, psi(1) = 0 | MATCH (relative difference at most 7e-16) |
| p = (delta/2) x wealth; bet = psi / sqrt((n - i + 1) v_(i-1)) (eq 1.4) | MATCH |
| Bet caps 1/m and 1/(1 - m), c = 1 | MATCH |
| Wealth freeze at 0 or 2/delta; strict acceptance below 2/delta on both tails (eq 1.6) | MATCH |
| Two one-sided tests intersected, total noncoverage delta (Thm 5.1) | MATCH |
| Bisection inversion (Prop 6.1: acceptance set is an interval) | MATCH in logic; 0 monotonicity violations in 200 random samples |
| Affine map to [-790, 790] | MATCH (same result as on the unit scale) |
| Early-return branch for m = 0 or 1 | Dead code (never fires); harmless |
| Midpoint acceptance check | Added safeguard, not in the paper |

Reviewer numeric checks: all-zeros n = 10 gives upper bound 0.3548, equal to the exact one-sided Clopper-Pearson bound. Coverage on iid Bernoulli data, 800 repeats per case: 0.9988, 0.975, 0.975 at nominal 0.95. No empty intervals.

**Reviewer verdict: PASS WITH NOTES.**

## Notes and the project's answer

| # | Reviewer note | Answer |
|---|---|---|
| 1 | High if it applies: the plain method is not valid for sampling without replacement from a fixed finite cohort (paper Sec 7 needs a different construction). | Does not apply. The interval's unit is the world. Each world is an independent draw from the declared world law with its own reserved seeds (`src/trackc_research_driver.py`, `build_schedule`), not a draw without replacement from a fixed list. Patients are averaged inside a world first. |
| 2 | Medium: validity needs n fixed before data and a fixed data order. | The amendment fixes n before any evaluation world runs. The analysis must feed world differences in world-index order (schedule order), never sorted or filtered. The analysis code must enforce this; a test is required. |
| 3 | Low: the +/-790 support had no source the reviewer found. | Source: utility is capped at 730 days and at worst loses 5 days x 12 courses = 60, so a paired difference lies in [-790, 790] (`tests/test_trackc_config_parity_extended.py`, `test_paired_difference_support_follows_from_utility_settings`). With the 10-day penalty sensitivity the lowest value is -120, still inside -790. |
| 4 | Low: dead branch. | Left in place; harmless. |
| 5 | Low: `None` (empty set) could be ignored silently. | The analysis must treat `None` as an inference failure and stop, never drop it. |

## Not verified

The reviewer did not re-derive the proofs (Thm 5.1, Prop 6.1, Prop 6.2). No human statistician has reviewed the method.
