# Track A TMZ development falsification 3 (frozen before run)

2026-11-23. Generated data only. Code `src/trackA_dev3_tmz.py` (reuses `src/trackA_dev2.py` arms unchanged). Output `output/trackA_dev3_tmz/result.json`. Seed 202611233. 20 replicates per scenario. 40 training origins, 100 evaluation origins.

## Why a redesign

Run 2 ([result](trackA_dev2_falsification_result_2026-11-23.md)) used fixtures where between-origin and within-origin effects were the same. Then any direction estimate finds the right genes, so shuffled links tie true links. True links can add accuracy only when an origin-level factor confounds the pooled fit and does not carry to new origins.

## Fixture

Each origin gets a hidden origin-level factor z. In training, z shifts gene 4 for all cultures of that origin and adds `confound * z` to their viability. Within an origin gene 4 noise is unrelated to viability. So true within-origin centering removes z. Pooled ridge and shuffled-link centering do not.

| Scenario | Training confound | Evaluation confound | Expected if mechanism is real |
|---|---|---|---|
| confounded_transfer | 0.8 | 0 | candidate beats base and shuffled links |
| confounded_both | 0.8 | 0.8 | base may win; candidate must not lose by > 0.01 |
| strong_confound | 1.5 | 0 | as confounded_transfer, larger gap |
| plain_transfer | 0 | 0 | run-2 control: tie with shuffled expected |
| null | 0 (no W effect) | 0 | no harm |

Within-origin effect W = (0.9, −0.6, 0, …) in all but null; between-origin effect = W.

## Arms and rule

Same arms as run 2: origin-balanced ridge (base), replicated direction (candidate), unchecked direction, shuffled-link direction, outcome-blind direction, one-slope contrast. Primary: per-origin absolute error of mean viability.

Pass (development only): in confounded_transfer AND strong_confound, replicate 95% interval of (base − candidate) > 0 AND of (shuffled − candidate) > 0. In confounded_both and null, mean (base − candidate) ≥ −0.01. Otherwise the sister-link gain claim stays unproven. This is a software and mechanism check, not evidence for biology.
