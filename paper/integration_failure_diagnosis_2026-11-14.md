# Integration: why the real-data checks failed (diagnosis, 2026-11-14)

Exploratory diagnosis. No locked test was re-run or changed.

## Code check (sign and feature drift)
Script: scratch only (frozen model refit on TCGA, same as confirm/UPENN scripts).
- TCGA in-sample Spearman rho of frozen pls2 MES model: 0.609.
- Dev repeated-CV rho (best of 300 tries): 0.354 (paper/integration_methoddev_lock_2026-11-12.md).
- Locked external rho: 0.114 (output/integration_methoddev/confirmation.json).
- No missing feature columns in CPTAC, IvyGAP or UPENN; NaN fraction <= 0.008; median per-feature shift 0.21-0.39 TCGA SDs.
- Reading: the drop 0.61 -> 0.35 -> 0.11 is the pattern of overfit plus selection (winner's curse), not a code bug. Feature drift is modest. Not fully ruled out: scanner harmonization (no ComBat used).

## Pattern across all four real-data MRI->biology tests
1. BAMD Gate 1, TCGA: failed (paper/integration_bamd_gate1_result_2026-11-12.md).
2. BAMD Gate 1b, Hu 2023 spatially matched biopsies: failed (paper/integration_bamd_gate1b_result_2026-11-12.md). This removes the "sampling mismatch" excuse.
3. Locked MES confirmation, CPTAC+IvyGAP: failed (confirmation.json).
4. UPENN survival: failed, opposite direction, post hoc (output/integration_methoddev/upenn_survival.json).
Reading: in these public data, routine MRI carries little transcriptomic-state signal. This is a method-premise problem: methods that need MRI to PREDICT biology will keep failing here.

## Prior art (haiku literature search, partly unverified)
Hu 2023 Nat Commun (spatial MRI-RNA, single site); Chaddad 2019 (radiomics+RNA AUC 0.779 vs RNA 0.731, unverified details); Kickingereder 2016/2018; Rathore 2018; ComBat leakage warning (Carre 2022). Gap seen: few studies test combined vs MRI-only vs RNA-only on the same score; no fusion that weights each source by its measured reliability in a treatment simulator was found (not a full prior-art search).
