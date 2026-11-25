# Track A TMZ: frozen repair-axis score on GSE7696 — RESULT (2026-11-24)

One locked run. Plan `trackA_tmz_gse7696_prereg_LOCKED_2026-11-24.md` (commit ee7f5c1). Code SHA-256 a6dc38d7...6584b3, unchanged since the lock commit. Output `output/trackA_tmz_gse7696/result.json`.

## Frozen rule: FAIL

68 primary GBM patients (TMZ+RT 42, RT 26), 58 deaths. Cox: OS ~ arm + score + arm x score + MGMT status + age.

| Quantity | Value |
|---|---:|
| arm x score coefficient | +0.403 (wrong direction; negative was predicted) |
| z | +1.53 |
| one-sided p | 0.94 |
| Ablation (MGMT expression alone) interaction z | +0.97 |
| Random 4-gene swaps: median z / 5th percentile | +0.45 / -0.96 |
| Share of random swaps with weaker z | 0.082 |

The score did not identify patients who benefit more from adding TMZ. The point estimate is in the opposite direction but not significant (z = 1.53).

## What this means

- The Mayo PDX pass (p = 0.040, n = 27) did **not** replicate in humans. Combined, the evidence does not support the repair-axis score as a TMZ-benefit predictor.
- The Mayo result should be reported as a single small-cohort preclinical finding that failed its independent human test.

## Weakest points

- Low power: 68 patients, 26 controls; an interaction test needs many more. A fail here means "not shown", but the direction was also wrong.
- Randomization of every GSE7696 sample is not verified.
- Different endpoint (overall survival in patients, mixing prognosis) and platform (microarray) from Mayo.
- PMS2 probes also match a pseudogene.

## Status of the Track A drug part

Positive: failed on the independent test. Unique: not established. Substantial: failed the shift test. GSE7696 and Mayo are now spent for this score. No retuning on either.
