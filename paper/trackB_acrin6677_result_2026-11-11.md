# ACRIN enhancing area forecast (2026-11-11)

## What was tested

The [amendment](trackB_acrin6677_amendment1_2026-11-11.md) switched the primary endpoint to central-reader enhancing area before reading measurements because the volume unit could not be established. The planning-time defect was that the history supplied one increment for two parameters. The amendment fixed growth to zero for falling histories and the decaying term to zero otherwise. The exploratory comparison used the first eligible window per patient. The [synthetic test record](trackB_acrin6677_test_record_2026-11-11.md) documents the checks before the single patient run.

## What was found

**Failure under the frozen rule.** In the [verified aggregate result](../output/track_b_acrin6677_2026-11-11/result.json), 85 patients were eligible. Median absolute log errors and percentile bootstrap 95% intervals were:

| Method | Median | Interval |
|---|---:|---:|
| B1 last value | 0.292 | 0.192–0.363 |
| B2 log-linear | 0.703 | 0.580–1.000 |
| B3 damped log-linear | 0.454 | 0.340–0.617 |
| M | 0.466 | 0.289–0.608 |

B1 was best. M minus B1 had median paired error difference +0.130, two-sided exact signed-rank p = 0.0000147, and M won 25 comparisons. Source for all values is the aggregate result above. The runner breaks equal baseline medians in protocol listing order and calculates the exact conditional sign distribution with zero differences dropped and tied absolute differences assigned average ranks. The [runner](../src/track_b_acrin6677_run.py) uses 10,000 patient bootstrap resamples with seed 20261112.

## Weakest point

Bevacizumab can reduce enhancement through changes in vascular permeability, so the MRI area can show pseudo-response. Each patient contributed a single window, and the baseline was selected using the same cohort. The nominal signed-rank p-value describes the comparison with that selected baseline. [Frozen protocol](trackB_acrin6677_protocol_2026-11-12.md)

## What comes next

Novelty was unchecked at the run. The [bounded prior-art search](trackB_acrin6677_prior_art_2026-11-11.md) records related work and remaining gaps. The [advanced-MRI draft](trackB_acrin6677_advanced_mri_draft_2026-11-11.md) is exploratory and unrun.
