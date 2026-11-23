# Track A TMZ development falsification 6 (frozen before run)

2026-11-23. Generated data only. Code `src/trackA_dev6_tmz.py` (fixtures and arms from dev4/dev5 unchanged). Output `output/trackA_dev6_tmz/result.json`. Fresh seed 202611236. 20 replicates. Same 7 scenarios.

## Candidate `halves_selector`

Inside training only: split origins into two fixed random halves. For each of {within_only, origin_balanced_ridge, mundlak_ridge}, fit on one half and score origin MAE on the other, both directions; sum the two. Pick the lowest sum. Refit the picked method on all training data.

## Known limit stated before the run

In confounded_transfer and confounded_both the TRAINING data have the same distribution (confound 0.8). Only the evaluation origins differ. A selector that sees only training data cannot tell these apart, so it cannot be best in both. This run measures how close it gets.

## Pass rule (development only)

Regret = candidate mean MAE − best single comparator mean MAE in that scenario (comparators: within_only, origin_balanced_ridge, mundlak_ridge, shuffled-link contrast). Pass if the candidate's maximum regret over all 7 scenarios is smaller than the maximum regret of every single fixed comparator. Otherwise it adds nothing over choosing one established method.
