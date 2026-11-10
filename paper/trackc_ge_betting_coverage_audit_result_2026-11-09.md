# Track C GE-betting audit and outcome-free coverage screen: result (2026-11-09)

**Status: synthetic methods result.** No Track C policy outcome was generated, no patient data were read, and no evaluation seed was reserved. The protocol ([file](trackc_ge_betting_coverage_audit_protocol_2026-11-09.md)) was committed in `f093344` before the screen ran. Result file: `output/trackc_ge_betting_coverage_audit_20261109/coverage.json`, produced in a scratch copy and copied unedited.

## 1. Audit of `src/trackc_ge_betting.py` against the published equations

Source: Martinez-Taboada and Ramdas, arXiv:2608.21694v2, read from the arXiv HTML on 2026-11-09 through an automated page query. The equations below were not re-typed by a human from the PDF; treat equation numbers as reported by that query.

| Paper item | Paper definition | Code | Verdict |
|---|---|---|---|
| Mean and variance path (Eq. 1.3) | `mu_t = (1/2 + sum X_j)/(t+1)`, `v_t = (1/4 + sum (X_j - mu_(j-1))^2)/(t+1)` | `_variance_path`, lines 38-52 | Matches |
| `psi(p) = phi(Phi^-1(p))/p` | Section 1 | `_psi`, lines 55-61 | Matches |
| Bet (Eq. 1.4) | `min{ psi((delta/2) K_(i-1)) / sqrt((n-i+1) v_(i-1)), c/m }`, lower tail `c/(1-m)`, default `c = 1` | `_tail_wealth`, lines 79-89 (0-based index, so `n - index = n-i+1`) | Matches |
| Absorption | wealth left unchanged at 0 or `2/delta` | lines 77-78, 91 | Matches; the code caps wealth at `2/delta` instead of keeping an overshoot, which gives the same accept/reject decision |
| Two-sided set (Eq. 1.6, `U+ = U- = 1`) | terminal `K+ < 2/delta` and `K- < 2/delta` | `_accepts_mean`, lines 95-100 | Matches |
| External randomization | optional | not implemented | Correctly excluded; it is not adopted |
| Interval shape (Theorem 5.1) | accepted set is an interval | endpoint bisection, lines 137-173 | The theorem covers the intersection. The code bisects each tail separately, which also assumes each one-sided set is monotone. New test `test_bisection_interval_matches_dense_grid_acceptance` compares the endpoints with a 401-point grid on 24 random samples and found no disagreement; the screen found 0 inversion inconsistencies in 12,000 intervals. This is empirical support, not a proof. |

**Assumptions of the validity theorem (Theorem 5.1):** observations in `[0, 1]` with `E[X_i | F_(i-1)] = mu`. Independence is not required. So validity for Track C rests on two project facts, not on the distribution shape:

1. **Support.** Every world-level paired difference lies in `[-790, 790]` days. This follows from the utility settings (horizon 730 days, 5 days per course, at most 12 courses); the new parity test `test_paired_difference_support_follows_from_utility_settings` pins the arithmetic. A course-penalty sensitivity of 10 days changes the support to `[-850, 850]` and needs its own interval support.
2. **Common mean across worlds.** Worlds must come from one fixed law, independently or at least with a common conditional mean. The seed schedule is still unmaterialized, so this is not yet demonstrated.

## 2. Screen result (30 cells, 2,000 repetitions each, same pass rule as the earlier screens)

| Method | Cells passed | Joint coverage range | Screen |
|---|---:|---:|---|
| GE-betting (deterministic, `c = 1`) | 30/30 | 95.95% to 100% | Passed |
| Hoeffding | 30/30 | 99.15% to 100% | Passed |
| Empirical Bernstein (Maurer-Pontil) | 30/30 | 100% | Passed |
| Student-t (reference) | 18/30 | 74.6% to 97.7% | Failed, as in the earlier screens |

GE-betting had 0 empty intervals and 0 inversion inconsistencies. Runtime was 14.3 seconds on 6 processes.

**Weakest point:** the pass uses six synthetic distribution pairs. It does not test the real world-difference law, which cannot be observed without policy outcomes. The GE-betting margin is thin in the maximum-variance cell (lowest Wilson lower bound 0.9499 at 100 worlds, two-point pair), which is consistent with an exact method at its limit.

## 3. Precision: the 14-day target was not met in any cell

Median of the larger family half-width, in simulated days:

| Worlds per family | GE-betting, SD 30 days | GE-betting, rare event | GE-betting, uniform | GE-betting, two-point | Hoeffding | Empirical Bernstein (SD 30) |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 400.5 | 584.0 | 667.5 | 790.0 | 724.1 | 790.0 |
| 14 | 224.1 | 353.1 | 410.1 | 601.4 | 625.1 | 790.0 |
| 30 | 113.3 | 176.4 | 243.7 | 374.9 | 427.0 | 663.7 |
| 100 | 36.8 | 60.2 | 117.9 | 189.1 | 233.9 | 198.9 |
| 207 | 18.3 | 31.6 | 78.0 | 128.6 | 162.6 | 97.6 |

The fraction of repetitions with both half-widths at most 14 days was 0 for every method in every cell. GE-betting is the narrowest valid method tested, by a wide margin. Even so, with a world-difference SD of only 30 days it reached 18.3 days at 207 worlds per family.

**Consequences, stated without changing any rule:**

- The earlier "at least 207 worlds per family" figure is only where one worst-case lower bound drops below 14 days. This screen shows that 207 worlds do not reach 14 days even in the most favourable cell tested.
- Extrapolating the SD-30 cell by `1/sqrt(n)` gives roughly 350 worlds per family for a 14-day median half-width. This is an extrapolation, not a tested cell, and it holds only if the real SD is about 30 days. The real SD is unknown; the required count scales with SD squared (about 1,400 worlds at SD 60, about 3,900 at SD 100).
- The current seed block allows 7 worlds per family. At 7 worlds the narrowest GE-betting median was 400.5 days.

## 4. Runtime estimate and its uncertainty

Only one timing observation exists: 1,373.3 s (single-response) and 1,482.7 s (mixture-response) for one world per family with two patients per cohort ([timing record](trackc_robust_method_research_and_timing_protocol_2026-11-06.md)). That is 0.793 hours per world pair.

| Worlds per family | Serial hours at 2 patients per cohort |
|---:|---:|
| 7 | 5.6 |
| 100 | 79.3 |
| 207 | 164.2 |
| 350 (extrapolated SD-30 need) | 277.7 |

Uncertainty: with one observation per family there is no measured run-to-run variance, so no statistical interval can be given. The two families differ by 8%. Cohort sizes are not yet locked; evaluation cost grows roughly in proportion to patients per cohort, so a cohort of 18 would multiply these figures by about 9. Treat every figure as a point estimate with unknown error, probably within a factor of two for the two-patient case and much larger once cohort size is set.

## 5. Decision

GE-betting is the strongest candidate found and passed its predeclared synthetic screen. It is **not adopted**: it still lacks independent statistical review and a materialized, independent world schedule. More importantly, no valid method tested can meet the unchanged 14-day precision target at any feasible world count on this machine. Confirmatory Track C evaluation therefore stays stopped. The honest report is a feasibility result: a valid interval method exists, and the precision requirement is not attainable with the available worlds and runtime.

A path that keeps every rule, for the user to decide: declare a fixed world count from a runtime cap before any outcome, use GE-betting at `delta = 0.025` per family, and report "inconclusive" when either half-width exceeds 14 days, as the config's `nonpass_rule` already requires. On present evidence that run would very probably end inconclusive, so it is not recommended unless the runtime budget grows to hundreds of hours.
