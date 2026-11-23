# Track A TMZ development falsification 5 (frozen before run)

2026-11-23. Generated data only. Adds the two comparators required by the [prior-art check](trackA_tmz_gated_contrast_prior_art_2026-11-23.md). Code `src/trackA_dev5_tmz.py` (fixtures and candidate unchanged from `src/trackA_dev4_tmz.py`). Output `output/trackA_dev5_tmz/result.json`. Fresh seed 202611235. 20 replicates. Same 7 scenarios as dev4.

## New comparators

- **mundlak_ridge**: origin-balanced ridge on [culture expression, origin-mean expression]. At prediction, the origin mean uses only the evaluation origin's own expression (no outcomes).
- **within_only** (analog of the family-GWAS direct-effect predictor): coefficients from ridge on within-origin centered expression and viability (repeated origins only); then a scalar slope and intercept fitted on training origin means (origin-balanced least squares) of x·β_within.

Kept: origin-balanced ridge, gated contrast (candidate), unchecked contrast, shuffled-link contrast.

## Pass rule (development only)

The strongest comparator in each confounded-transfer scenario is the best of {origin-balanced ridge, mundlak_ridge, within_only, shuffled-link contrast}. In confounded_transfer, strong_confound and unequal_confound: replicate 95% interval of (strongest comparator − candidate) > 0. In reversal, confounded_both and null: candidate mean MAE ≤ best comparator mean + 0.01. Any miss means the candidate adds nothing beyond established within/between methods in these fixtures.
