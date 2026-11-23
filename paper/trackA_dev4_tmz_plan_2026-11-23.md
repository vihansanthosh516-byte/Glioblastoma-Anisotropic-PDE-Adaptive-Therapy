# Track A TMZ development falsification 4 (frozen before run)

2026-11-23. Generated data only. Code `src/trackA_dev4_tmz.py`. Output `output/trackA_dev4_tmz/result.json`. Seed 202611234 (fresh). 20 replicates. 40 training origins, 100 evaluation origins.

## Candidate `gated_contrast`

- Model: existing one-slope contrast (`_fit_shared`, within-origin weight 1, ridge penalty 0.2).
- Gate inside training only: split training origins into two fixed random halves. Fit contrast and origin-balanced ridge on half A, score origin MAE on half B; then the mirror. Admit the contrast if the SUM of the two improvements (ridge MAE − contrast MAE) is > 0. Otherwise use origin-balanced ridge.
- The sum rule is looser than the run-3 "both halves" rule, which admitted only ~50% of signal draws. This change was chosen after seeing run 3 and is labelled post hoc.

## Arms

origin-balanced ridge (base); gated_contrast (candidate); contrast unchecked; shuffled-link contrast (links shuffled in the contrast term only, whole rows kept); replicated_direction (run-2/3 candidate).

## Scenarios

| Scenario | Between-origin effect | Train confound | Eval confound | Sisters |
|---|---|---|---|---|
| reversal | −W | 0 | 0 | 1–4 |
| confounded_transfer | W | 0.8 | 0 | 1–4 |
| strong_confound | W | 1.5 | 0 | 1–4 |
| confounded_both | W | 0.8 | 0.8 | 1–4 |
| plain_transfer | W | 0 | 0 | 1–4 |
| unequal_confound | W | 0.8 | 0 | 1–8 |
| null | 0 (and within 0) | 0 | 0 | 1–4 |

## Pass rule (development only)

In confounded_transfer, strong_confound and unequal_confound: replicate 95% interval of (base − candidate) > 0 and of (shuffled − candidate) > 0. In reversal, confounded_both and null: mean (base − candidate) ≥ −0.01. Any miss = fail. Software and mechanism check only, not biology.
