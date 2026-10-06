# Review packet: GBM forecasting project, state on 2026-10-05 (Phases 0-2 done, Phase 3 started)

Please check this critically. Every number cites a file in the repo. Interim items are labelled INTERIM. Negative findings are included on purpose.

## 1. What the project is
A computational study of glioblastoma (GBM) tumor-change forecasting. Question: does a PDE growth model (Fisher-KPP reaction-diffusion, with or without white-matter tract anisotropy and treatment terms) forecast the next MRI tumor mask better than simple baselines, on patients it was not fitted on?
Data: MU-Glioma-Post (203 patients; 134 eligible; MRI masks, 2 mm grid, target = tumor core). External test set: LUMIERE (91 GBM patients; only volume tables and RANO ratings so far; the zip with images and automated masks was downloaded and its MD5 checked, not yet unzipped).
Everything is simulation or retrospective. Nothing is clinical.

## 2. What was done (Phases 0-2, plus start of 3)
- **Preregistration-style plan** `analysis_plan_v1.md` with Amendments 1-3 (written before the results they govern). Primary endpoint: per-patient mean Dice gain vs persistence (copy the last mask forward). "Broadly positive" needs BOTH mean CI lower bound > 0 AND median > 0. Holm / BH-FDR for multiple tests.
- **Manifests and leakage tests** (`data/manifests/`, `tests/test_no_leakage.py`): patient-level 5-fold split by hash; 203 MU patients, 134 eligible (104 GBM, 30 other), 334 primary pairs (`consort_mu.json`).
- **Track A audit** (`output/trackA_audit.json`, `output/cgat_leak_free.json`): single-cell zone classifier.
- **Cohort tables** (`output/cohort_tables.json`).
- **Solver verification** (`output/solver_verification.json`).
- **Baseline ladder** (`output/baseline_ladder.json`), **forecastability** (`output/forecastability.json`), **DTI alignment** (`output/dti_alignment.json`), **subgroups** (`output/subgroups.json`).
- **Repo repair**: 5 old failing tests fixed; 93 tests pass.
- **PDE re-run on the manifest split**: INTERIM, still computing.

## 3. Results, each with its weakest point

**3.1 Track A (cell zone classifier) is weak** (ledger 2o; `trackA_audit.json`, `cgat_leak_free.json`).
- Leak-free graph classifier: 63.8% accuracy (SD 13.3) over 5 patient-level folds. A model that only sees cells-per-patient gets 64.6% (SD 27.5). PCA32+logistic regression gets 66.5%.
- The "Healthy" class is 3 donors with other diagnoses (neurocytoma, oligodendroglioma, meningioma). 13 of 21 patients contribute one class only.
- Earlier higher accuracies (about 75%) leaked through label edges and cVAE. Track A is exploratory only.
- Weakest point: 21 patients.

**3.2 MU cohort selection is informative** (ledger 2p; `consort_mu.json`).
- Eligible patients progressed 88.8% vs 47.8% of excluded (SMD 0.97). The forecastable set is mostly progressors.
- MU vs LUMIERE shift: median age 58 vs 61; MGMT methylated 42% vs 56%; median scan interval 77 vs 91 d.
- Weakest point: results may not transfer to stable patients.

**3.3 Solver passes numerical checks** (ledger 2q; `solver_verification.json`; synthetic).
- Spatial order about 2.0; temporal order about 1.0-1.1; mass change 1e-8; domain boundary respected.
- Found: with strong anisotropy the u>=0 clamp adds 2.8% mass in a stress test. Fixed: clamp mass is now logged.
- Found: at the 2 mm production grid the front speed is off by +50% at (D 0.01, rho 0.1) vs a 0.25 mm reference; 7 of 12 cells within 10%. At the 70-day horizon the radius error is at most 0.64 mm (pre-declared limit 1 mm).
- One test (V7a) was invalid; V7b was added after the first run and before its results were read. This is disclosed.
- Weakest point: verification is not validation of the biology; the 1 mm limit is loose next to the Dice effects.

**3.4 Baseline ladder** (ledger 2r; `baseline_ladder.json`; 133 patients, 328 pairs; INTERIM 2,000-draw bootstrap, the primary run will use 10,000).
- Persistence mean Dice 0.277.
- Delta vs persistence (CI): geometric growth rule, no PDE, +0.0136 (+0.0096 to +0.0173), better in 68%, broadly positive; last_rate -0.0445; linear -0.028; gompertz -0.0026 (CI includes 0).
- Oracle that uses the true next volume (upper bound, not a baseline): +0.128.
- **CORRECTION:** an earlier ledger entry said +0.031 for the geometric rule. That did not match the file. I fixed the ledger, plan A3.2 and W23. Please check the ledger yourself.
- Weakest point: the geometric rule uses brain mask and distance ordering, not biology.

**3.5 Forecastability is poor** (ledger 2s; `forecastability.json`).
- 57% of pairs have change below the assumed noise floor (SNR < 1). A shift of one voxel at the boundary gets Dice 0.53; persistence gets 0.28.
- Weakest point: the noise floor is an assumption and probably too high at 2 mm.

**3.6 DTI direction test is null** (ledger 2t; `dti_alignment.json`; exploratory).
- Growth direction vs atlas principal direction: |cos| 0.532 (CI 0.494-0.569). Permutation null mean 0.515; p = 0.32. 124 patients.
- Weakest point: a population atlas, not patient DTI; centroid shift is a coarse direction. This does not prove tracts are irrelevant.

**3.7 Subgroups, bad misses, residuals** (ledger 2v; `subgroups.json`; exploratory, baselines only, INTERIM bootstrap).
- Geometric-rule gain is +0.010 to +0.025 in every stratum with enough patients (prior growth, days since radiotherapy, horizon). "Stable" has 10 patients and its CI includes 0.
- Bad-miss rate (Dice < 0.1 when persistence >= 0.1): persistence 0%, geometric 0.4%, last_rate 13.9%, linear 16.5%, gompertz 2.6%.
- Gain is not related to volume, interval, prior growth or days since radiotherapy (all |rho| <= 0.06).
- Weakest point: the gain is tiny; the 0.1 threshold is analysis-defined, not clinical.

**3.8 Tests**: 93 pass (`tests/`). 5 old failures were fixed: two code restores, one API fix, one assertion rewritten to record that D is not identifiable (D RMSE 93% relative).

## 4. INTERIM / not done
- PDE re-run (script 100): two-stage design (full 125-cell grid on each patient's first pair, then chosen cell per fold on later pairs). About 7 of 127 first pairs cached when I stopped; about 16 h expected. No PDE result exists yet. Key test when done: PDE vs geometric rule, paired per patient.
- Not run yet: 1 mm grid pilot, compartment forecasts, negative controls, segmentation perturbation, freeze (`configs/development.yaml`), one external run on LUMIERE, failure atlas, 10,000-draw primary run.
- LUMIERE masks are machine-made (HD-GLIO-AUTO, DeepBraTumIA), not expert masks.

## 5. Open weak points (`vault/WEAK_POINTS.md`, W1-W26)
Main ones: cohort has non-GBM patients (GBM-only is a pre-set sensitivity run); treatment constants alpha 0.08 (TMZ) and beta 0.03 (RT) are assumed; the grid-edge fit; one 2 mm resolution; Track A; interim bootstraps; the user still has to fill `docs/isef_compliance.md`.

## 6. Questions for the reviewer
1. Is "better than persistence" a fair primary endpoint, given the oracle and geometric rule?
2. Is the two-stage grid selection leak-free (cells chosen on training patients' first pairs only)?
3. Is a +0.0136 Dice gain worth any claim?
4. What would you test next?
Also: please spot-check ledger numbers against `output/*.json`. One error was already found this way.
