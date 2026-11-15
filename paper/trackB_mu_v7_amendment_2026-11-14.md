# Track B v7 (MRI brightness + composition) — amendment, locked before v7 scoring (2026-11-14)

## Unchanged
Same question, data, pairs, folds, target, primary score, 6 baselines, bootstrap seed (`paper/trackB_mu_tcks_protocol_2026-11-14.md`). Multiple-try rule: try 7 passes only if the upper bound of the 99.29% (1 - 0.05/7) CI vs best baseline is < 0.

## v7 method (only change)
All v5 features plus brightness from the INPUT scan only (t1c and t2f images at tp_in). Each image is divided by the median of its nonzero (brain) voxels:
- t1c_enh: mean t1c in enhancing label 3;
- t1c_core_sd: std of t1c in core (labels 1+3);
- flair_edema: mean t2f in edema label 2;
- flair_shell: mean t2f in a 5-voxel shell just outside the whole tumor (dilation by 5, minus tumor, inside brain).
Empty region -> 0. Same model as v5: median (L1) linear regression, standardized on training folds. No tuning. v6 edge-shape features are not used.

## Honesty note
Post hoc, try 7, reused cohort. Exploratory only. Intensity radiomics for glioma growth has prior art.
