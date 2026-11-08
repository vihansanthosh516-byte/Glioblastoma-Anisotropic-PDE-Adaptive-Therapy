# Track C max-T coverage audit result

**Status:** predeclared outcome-blind synthetic screen failed. The tested observed-SE max-T interval method is not eligible for research-scale inference or a positive claim.

## Result

The audit used six fixed distribution pairs, four world counts (5, 10, 14, 20), 2,000 Monte Carlo repetitions per cell, and 2,000 bootstrap resamples per repetition. All 24 cells completed with zero inference failures. The predeclared rule required every joint-coverage Wilson 95% lower bound to be at least 0.93. No cell met that rule.

| Worlds per family | Mean joint coverage across six pairs | Lowest cell coverage | Lowest Wilson lower bound |
| ---: | ---: | ---: | ---: |
| 5 | 75.0% | 70.05% | 68.01% |
| 10 | 86.0% | 82.75% | 81.03% |
| 14 | 89.0% | 86.15% | 84.58% |
| 20 | 91.0% | 88.7% | 87.22% |

Lowest coverage occurred for centered Gamma(shape 2) at 5 worlds: **70.05%** (Wilson 95% interval 68.01%–72.02%). Highest coverage was **93.2%** for Student t(df 3) at 20 worlds (Wilson 95% interval 92.0%–94.2%). These are below the nominal 95% target. The two-lower-bounds-above-zero rule produced 41 false positives in 48,000 global-null repetitions (0.085%); this does not offset undercoverage.

## Interpretation and action

The helper's intervals are too narrow under the tested conditions, especially with few worlds and skewed outcomes. Do not use this helper for a Track C policy evaluation, a positive decision, or a claim of superiority. Do not tune it using archived or future policy outcomes.

Next, obtain independent statistical review and predeclare a replacement interval or test, including its estimand, dependence structure, precision target, and failure rule. Validate the replacement with fresh outcome-blind synthetic draws matched to the world-level data structure. If no method meets its coverage and precision criteria at feasible world counts, report descriptive simulation estimates with uncertainty limits stated plainly, or stop confirmatory comparison. Keep research-scale evaluation blocked until the replacement method passes.

## Provenance

- Fixed protocol: [`trackc_max_t_coverage_audit_protocol_2026-11-08.md`](trackc_max_t_coverage_audit_protocol_2026-11-08.md).
- Implementation: `src/trackc_max_t_coverage_audit.py`.
- Machine-readable output: [`../output/trackc_max_t_coverage_audit_20261108/coverage.json`](../output/trackc_max_t_coverage_audit_20261108/coverage.json).
- Runtime: 22.32 seconds; Python 3.14.4; NumPy 2.4.4.
- The output's code, helper, protocol, and integrated-config fingerprints matched the files present when the audit ran. They were checked after the run and before the config status update.
- No patient data or Track C policy outcomes were used. No figure changed.
