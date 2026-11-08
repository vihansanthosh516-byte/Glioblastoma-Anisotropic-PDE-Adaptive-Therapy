# Track C Student-t Bonferroni candidate screen result

**Status:** candidate failed its predeclared outcome-blind screen. Do not use it for policy evaluation or a positive claim. Track C inference remains blocked.

## Result

The audit tested two-sided Student-t intervals with a Bonferroni adjustment across the two primary response-family means. It used six fixed distribution pairs, 20,000 Monte Carlo repetitions per pair and world count, and world counts 5, 10, 14, and 20. All 24 cells completed with zero inference failures. No patient data or Track C policy outcomes were used.

| Worlds per family | Cells passing screen | Lowest joint coverage | Lowest Wilson 95% lower bound |
| ---: | ---: | ---: | ---: |
| 5 | 2/6 | 89.8% | 89.4% |
| 10 | 4/6 | 91.0% | 90.1% |
| 14 | 4/6 | 91.0% | 90.9% |
| 20 | 4/6 | 92.0% | 91.8% |

Across all cells, joint coverage ranged from **89.815%** (centered Gamma(shape 2) in both families, 5 worlds) to **96.505%** (standardized Student t(df 3) in both families, 5 worlds). The predeclared rule required every Wilson lower bound to be at least 93%; only 14 of 24 cells passed. The largest observed both-lower-bounds-above-zero rate was 7/20,000 in one cell; the total was 40/480,000. This low false-positive count does not remedy undercoverage.

## Interpretation

Bonferroni controls simultaneous coverage when each family-specific interval has its stated marginal coverage. The Student-t interval is exact only for independent normally distributed world differences within each family. The centered-Gamma results show that this candidate does not meet the screen under the tested skewed distributions at the proposed world counts. Do not treat it as a validated replacement for the failed observed-SE max-T method.

Primary utility lies in `[-60, 730]` simulated days, so a paired candidate-control world difference lies in `[-790, 790]`. A two-family, 95% Hoeffding bound has a half-width of about **625 days at 14 worlds**. A worst-case half-width of 14 days would require about **27,900 independent worlds per family**. This bound is conservative, but it shows that the distribution-free route cannot meet the current precision target at the existing 14-world seed capacity.

## Next gate

Do not run research-scale policy comparisons yet. Choose a separate, dated plan for either (1) more independent worlds, a runtime cap, and a fresh coverage/precision validation under the actual world-generation design, or (2) a descriptive simulation report that drops confirmatory coverage and superiority language. Do not select a method by looking for a positive policy result. Keep the finite-set minimax contrast and all Track A/Track B records unchanged.

## Provenance

- Protocol: [`trackc_t_bonferroni_coverage_protocol_2026-11-08.md`](trackc_t_bonferroni_coverage_protocol_2026-11-08.md).
- Implementation: `src/trackc_t_bonferroni_coverage_audit.py`.
- Machine-readable result: [`../output/trackc_t_bonferroni_coverage_audit_20261108/coverage.json`](../output/trackc_t_bonferroni_coverage_audit_20261108/coverage.json).
- Synthetic helper tests: 8 passed in a scratch copy using Python 3.14.4, NumPy 2.4.4, and SciPy 1.18.0.
- Code and protocol SHA-256 values in the output match the files used for the run.
- No policy outcomes, patient data, or sealed seeds were used.
