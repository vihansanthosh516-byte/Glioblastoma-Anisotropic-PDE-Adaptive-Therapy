# Integration method-dev lock (written 2026-11-13, before any CPTAC-GBM or IvyGAP feature is read)

Search finished: `output/integration_methoddev/log.json` (300 tries, TCGA n = 91 with RNA + MRI; sha256 1514dc4b...).
TCGA is spent. The best of 300 tries is optimistic (selection bias).

## Locked choice (one target, one method)

- Target: `neftel_MES` (Neftel mesenchymal-like state score from tumor RNA).
- Features: `new` (73 features from `src/integration_methoddev_features.py`).
- Model: `pls2` (PLS regression, 2 components), fit on all 91 TCGA patients, then frozen.
- Dev CV rho (10 x 5): 0.354 (sd 0.051). Optimistic; expect lower on new data.
- Direction for UPENN backup: higher score = higher hazard (target is MES).

## Confirmation rule (unchanged from `paper/integration_protocol_methoddev_2026-11-12.md`)

Apply the frozen model to CPTAC-GBM (23) and IvyGAP (26). Fisher-z combined Spearman rho >= 0.25 and p < 0.05. One run.

## Code and data hashes (sha256)

- src/integration_methoddev_features.py: bc6839d91087e8493fce0b749d58062989bfa220e00c7677d597bd9e9ce47515
- src/integration_methoddev_search.py: 49751ba50151acb52dc890dca0d71c8e29ffa1ae203cc03ed75161fedab7d4d8
- output/integration_methoddev/log.json: 1514dc4bdca80da46a8164cb5aed3cec4364c382e95b465eb97eb34c8eaa8bd2
- output/integration_methoddev/features_TCGA-GBM.csv: b07f51bbae2daa890e3c604508e43d2e075a4cdb401d2b543d525765b5398641

## Weak points

Small confirm set (49), so power is low. MES is the only target that held up; other targets were weak in dev.
