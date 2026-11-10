# Track C GE-betting outcome-free coverage screen: protocol (2026-11-09)

**Status: predeclared before the screen runs.** This protocol was written and committed before `src/trackc_ge_betting_coverage_audit.py` was executed. It uses synthetic bounded numbers only. It generates no Track C policy outcomes, reads no patient data, and reserves no evaluation seeds. Passing this screen does not adopt any interval method and does not authorize policy evaluation.

## Question

Under independent, identically distributed, bounded world-level differences on the declared support `[-790, 790]` simulated days, does each candidate fixed-horizon interval reach at least the familywise 95% coverage required for the two primary response-family contrasts, and how wide is it?

## What this screen can and cannot show

The GE-betting theorem (Martinez-Taboada and Ramdas, arXiv:2608.21694v2, Theorem 5.1, as read from the arXiv HTML on 2026-11-09) requires observations in `[0, 1]` with a common conditional mean, `E[X_i | F_(i-1)] = mu`; independence is sufficient but not required. Its validity therefore does not depend on the shape of the world-difference distribution. It depends on two project facts that this screen cannot test:

1. every world-level paired difference lies in `[-790, 790]` days (a property of the utility definition and the driver);
2. worlds are generated independently from one fixed law (a property of the materialized seed schedule and sampler).

The real world-difference distribution can only be observed by running policies, which stays blocked. This screen therefore checks the implementation and compares methods on declared synthetic distributions, including adversarial ones. It is not a validation under the project's actual world law.

## Methods compared

All intervals are two-sided for one family mean at `delta = 0.025`, so two families give familywise 95% by Bonferroni.

- `ge_betting`: deterministic GE-betting, cap `c = 1`, no external randomization (`src/trackc_ge_betting.py`).
- `hoeffding`: half-width `1580 * sqrt(ln(2/delta) / (2n))`.
- `empirical_bernstein`: Maurer-Pontil (2009) two-sided bound on the unit scale, `sqrt(2 V ln(4/delta) / n) + 7 ln(4/delta) / (3 (n - 1))`, with `V` the sample variance, clipped to the support.
- `student_t`: ordinary t interval. Reference only; it was already rejected by earlier screens and is not a candidate.

## Cells

World counts per family: `7, 14, 30, 100, 207`.

Distribution pairs (family 1, family 2), all on `[-790, 790]` days:

| Pair | Family 1 | Family 2 |
|---|---|---|
| `rare_pos_rare_pos` | `+790` with probability `p_n`, else 0 | same |
| `rare_pos_rare_neg` | as above | `-790` with probability `p_n`, else 0 |
| `two_point_two_point` | `+790` or `-790`, each 0.5 | same |
| `uniform_uniform` | uniform on the support | same |
| `narrow_normal_narrow_normal` | normal, mean 10, SD 30 | same |
| `skew_beta_narrow_normal` | `1580 * Beta(0.5, 5) - 790` | normal, mean 10, SD 30 |

`p_n = 1 - 0.025^(1/n)`, the rare-event probability used in the project's all-zero-sample lower bound. That gives 5 x 6 = 30 cells.

## Fixed settings

- Monte Carlo repetitions per cell: 2,000.
- Random streams: `PCG64DXSM(SeedSequence([20261109, cell_index, family_index]))`.
- Coverage of `ge_betting` is evaluated by testing acceptance of the true mean in every repetition. The full interval is computed on the first 200 repetitions of each cell for width, and on those repetitions the interval must contain the true mean exactly when the true mean is accepted (tolerance `1e-8` on the unit scale). Any disagreement fails the cell.

## Pass rule (same criterion as the earlier Track C screens)

A method passes a cell when the 95% Wilson lower bound of its joint coverage (both family intervals contain their true means) is at least 0.93, with no inference failure (empty interval). A method passes the screen only when it passes all 30 cells. This rule is not changed after results are seen.

## Reported precision (descriptive, no pass rule)

For each method and cell: median and 90th percentile of the larger of the two family half-widths, and the fraction of repetitions in which both half-widths are at most 14 days. The 14-day target is unchanged. These numbers describe synthetic distributions only; they do not predict the half-width of a real evaluation, because the real world-difference variance is unknown without policy outcomes.

## Outputs

`output/trackc_ge_betting_coverage_audit_20261109/coverage.json`, produced in a scratch copy and copied unedited into the repository.
