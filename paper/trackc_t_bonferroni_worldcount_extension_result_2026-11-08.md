# Track C Student-t Bonferroni world-count extension result

**Status:** the candidate failed the predeclared larger-count screen. It remains ineligible for policy evaluation. No patient data or Track C policy outcomes were used.

## Result

The interval method stayed fixed: two-sided Student-t intervals with Bonferroni adjustment across the two primary response families. The extension tested 11 distribution pairs at 30, 50, 100, and 200 worlds per family, with 20,000 Monte Carlo repetitions per cell. All 44 cells completed with zero inference failures.

| Worlds per family | Cells passing | Lowest joint coverage | Lowest Wilson 95% lower bound |
| ---: | ---: | ---: | ---: |
| 30 | 8/11 | 91.835% | 91.447% |
| 50 | 9/11 | 93.135% | 92.776% |
| 100 | 8/11 | 88.480% | 88.030% |
| 200 | 10/11 | 92.390% | 92.014% |

Across all cells, coverage ranged from **88.48%** to **95.9%**. The weakest cell was standardized Bernoulli(p = 0.1) in both families at 100 worlds. Coverage was not monotonic in world count: this cell fell well below the neighboring counts. The predeclared screen required every Wilson lower bound to be at least 93%; only 35 of 44 cells passed. The largest global-null rate for both lower bounds above zero was 5/20,000 in one cell. Low false-positive frequency does not compensate for undercoverage.

The additional bounded-skew distributions exposed a weakness that the initial six-pair screen did not cover. More worlds alone did not make this candidate pass the all-cells rule. Do not choose a world count from this grid or use t-Bonferroni for the integrated positive rule.

## Runtime feasibility

The prior timing-only pilot used two patients in each of the selection, calibration, and evaluation cohorts. It measured 1,373.3 seconds for one single-response world and 1,482.7 seconds for one mixture-response world. Linear extrapolation gives about **79.3 hours for 100 worlds per family** and **158.7 hours for 200 worlds per family**, even at those two-patient cohort sizes. These are rough projections from one timing observation per family. Larger patient cohorts would increase the work; exact scaling is unmeasured.

## Next gate

The observed-SE max-T and Student-t Bonferroni candidates both failed their fixed coverage screens. A distribution-free Hoeffding interval remains mathematically valid under independent bounded worlds but is far wider than the 14-day precision target at the current world capacity. Do not run research-scale policy comparisons. The next defensible choice is a new statistical design with independent review and feasible precision/runtime validation, or a descriptive simulation report that drops confirmatory superiority claims. Do not tune any method against policy outcomes.

## Provenance

- Protocol: [`trackc_t_bonferroni_worldcount_extension_protocol_2026-11-08.md`](trackc_t_bonferroni_worldcount_extension_protocol_2026-11-08.md).
- Implementation: `src/trackc_t_bonferroni_worldcount_extension.py`.
- Prior helper: `src/trackc_t_bonferroni_coverage_audit.py`.
- Machine-readable record: [`../output/trackc_t_bonferroni_worldcount_extension_20261108/coverage.json`](../output/trackc_t_bonferroni_worldcount_extension_20261108/coverage.json).
- Scratch-copy tests: 4 passed in 7.24 seconds using Python 3.14.4, NumPy 2.4.4, and SciPy 1.18.0.
- Code, helper, and protocol SHA-256 values match the files used for the run.
- No patient data, policy outcomes, or sealed seeds were used.
