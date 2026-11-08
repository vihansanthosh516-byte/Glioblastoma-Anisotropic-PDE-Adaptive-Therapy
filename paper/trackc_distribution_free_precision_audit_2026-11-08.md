# Track C distribution-free precision audit

**Status:** outcome-free mathematical and resource audit. It does not adopt an interval method or authorize policy evaluation. No patient data or policy outcomes were used.

## Bound from the current utility

The machine-readable protocol sets a 730-day progression horizon, 12-course maximum, and a five-day utility penalty per course. `robust_utility` takes the minimum of two progression times capped at 730 days, then subtracts the course penalty. `_checked_outcome` requires nonnegative times and course counts, and limits courses through the cutoff to total courses, which cannot exceed the configured budget.

Therefore each patient's utility lies in `[-60, 730]` days. A candidate-minus-control patient difference lies in `[-790, 790]`. An equally weighted mean of paired patient differences within one world has the same bounds. The range width for one independent world-level observation is therefore `R = 1,580` days.

This uses only the two currently named primary response-family contrasts. More confirmatory contrasts or anchor-specific families would require a stricter multiplicity allocation.

## Simultaneous Hoeffding intervals

For two response-family means and 95% familywise confidence, Bonferroni assigns 0.025 failure probability to each two-sided interval, or `delta = 0.0125` to each tail. Hoeffding's inequality gives the conservative half-width

`h(n) = R * sqrt(log(1 / delta) / (2 * n))`.

This finite-sample bound needs independent bounded world-level observations within each response family. It does not need normality. It does not require independence between the two response families because the two-sided familywise guarantee uses Bonferroni.

| Worlds per family | Hoeffding half-width |
| ---: | ---: |
| 7 | 883.96 days |
| 14 | 625.05 days |
| 30 | 426.99 days |
| 50 | 330.75 days |
| 100 | 233.87 days |
| 200 | 165.37 days |

The 14-day target requires **27,907 worlds per family** under this bound. That is a calculated threshold for this bound, not a proposed sample size and not a universal lower bound for every inferential method. The current development seed block has 89 free scenario seeds and supports at most 14 total worlds under the six-scenario-seeds-per-world schedule. No seeds were reserved.

## Runtime implication

The timing-only record measured 1,373.3 seconds for one single-response world and 1,482.7 seconds for one mixture-response world, each with only two patients in each of the selection, calibration, and evaluation cohorts. Linear extrapolation to 27,907 worlds per family is about **22,139 wall-clock hours (922 days)** for both families at those tiny cohort sizes. This is a rough extrapolation from one timing observation per family, not a bound; larger cohorts add work.

## Alternative checked

Maurer and Pontil's empirical-Bernstein theorem uses bounded observations and sample variance. Its additive range correction alone, rescaled to `R = 1,580` and the same per-tail `delta = 0.0125`, is 1,439.27 days at 14 worlds and 94.02 days at 200 worlds. It therefore cannot meet the 14-day target at either current sample-size reference, even if its variance term were zero. This rules out that particular bound for the current target at `n <= 200`; it does not rule out all variance-sensitive methods.

## Decision implication

The current 14-day precision target, broad utility range, seed capacity, and timing do not fit a distribution-free mean analysis. Do not substitute Student-t or observed-SE max-T: both failed their fixed outcome-blind screens. A tighter interval would need a defensible additional assumption or a different estimand, fixed before policy outcomes are inspected and reviewed independently. Otherwise, report the locked finite-world benchmark descriptively and make no population-level superiority claim. Do not choose a method because it produces a positive result.

## Sources and reproducibility

- Hoeffding, W. (1963). Probability inequalities for sums of bounded random variables. *Journal of the American Statistical Association, 58*(301), 13–30. [https://doi.org/10.1080/01621459.1963.10500830](https://doi.org/10.1080/01621459.1963.10500830)
- Maurer, A., & Pontil, M. (2009). Empirical Bernstein bounds and sample variance penalization. *Proceedings of the 22nd Annual Conference on Learning Theory*. [Primary paper](https://www.learningtheory.org/colt2009/papers/012.pdf)

Local inputs checked: `src/211_trackc_policy_benchmark.py` (`PolicyBenchmarkSettings`, `_checked_outcome`, `robust_utility`); `configs/trackABC_integrated_v2.json`; `src/trackc_seed_registry.py`; and `output/trackc_robust_policy_timing_pilot_20261106_followup/manifest.json`. The registry check read free-seed status only. It did not reserve seeds. Timing inputs were the saved timing record only.
