# Track B v6 (composition + edge shape) — amendment, locked before v6 scoring (2026-11-14)

## Unchanged
Same question, data, pairs, folds, target, primary score, 6 baselines, bootstrap seed (`paper/trackB_mu_tcks_protocol_2026-11-14.md`). Multiple-try rule: try 6 passes only if the upper bound of the 99.17% (1 - 0.05/6) CI vs best baseline is < 0.

## v6 method (only change)
All v5 features (input-scan mask only) plus edge shape of the input core (labels 1+3):
- core surface ratio: boundary voxels (core minus its 1-voxel 6-connected erosion) / core voxels;
- sphericity: pi^(1/3) (6 V)^(2/3) / A, with V = core voxels, A = boundary voxel count;
- edema reach: log((edema + core voxels) / core voxels).
Same model as v5: median (L1) linear regression, standardized on training folds. No tuning.

## Honesty note
Post hoc, try 6, reused cohort. Exploratory only. Shape radiomics for glioma growth has prior art.
