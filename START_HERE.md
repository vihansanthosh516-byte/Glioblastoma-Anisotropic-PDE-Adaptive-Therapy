# START HERE

Updated 2026-10-22. Read this file before continuing work. The main project is now the earlier GBM Tracks A, B, and C work for a local science fair.

## 1. Current goal

Begin coding the GBM Tracks A, B, and C project from the final proposed plan. Keep all three tracks equally important and separate. Do not treat the tracks as one fused patient-level model.

The current forward plan is `paper/research_plan_complete_project_form.md`, with details in `paper/research_plan_complete_project_appendix.md`. Tracks A, B, and C are equally important. The plan contains questions and procedures, not project outcomes. Each analysis remains conditional on its data, rights, and protocol gates.

## 2. Project shape

- **Track A — molecular and spatial data:** single-cell and spatial GBM analyses, including cell-state or tumor-zone prediction.
- **Track B — tumor growth:** longitudinal MRI, reaction-diffusion/PDE models, parameter estimation, and forecasting.
- **Track C — treatment suppression:** adaptive therapy, optimal control, and reinforcement-learning simulations.

These tracks share a GBM focus, but use different cohorts and endpoints. The repository does not establish a patient-level join from Track A to the MRI cohorts. Do not describe the existing work as one fused patient-level model or as a validated digital twin.

## 3. What the old evidence says

- **Track A:** the random-cell split leaked patient identity. The recorded patient-level C-GAT accuracy is 0.638, close to the cells-per-patient baseline of 0.646. Donor composition remains a confound.
- **Track B:** existing PDE forecasts did not demonstrate an advantage over simple baselines. The fitted diffusion parameter hit a bound in most cases in the recorded inverse problem.
- **Track C:** treatment suppression results are simulations using assumed drug-effect values. The “real_test” policy result uses fitted growth rates but assumed kill rates; it is not a real-drug evaluation.
- **Cross-track link:** the recorded Track A-to-B score was constant, and the source lineage is unresolved. No molecular-to-MRI patient linkage is established.

Use [README_STATUS.md](README_STATUS.md) for detailed numbers and the [findings ledger](vault/PAPER_FINDINGS_LEDGER.md) for cited result files. These are historical results, not proof that a new model works.

## 4. Working index

- [GBM project file guide](science_fair_gbm/README.md)
- [Track A/B/C improvement backlog](paper/GBM_TRACKS_IMPROVEMENT_BACKLOG.md)
- [Historical result status](README_STATUS.md)
- [Data-source inventory](vault/PAPER_DATA_SOURCES.md)
- [Old-track audit](paper/codex_audit/old_tracks.md)
- [Findings ledger](vault/PAPER_FINDINGS_LEDGER.md)
- [Archived GBM documents](archive_old_gbm/README.md)

Prior drafts exist. `paper/research_plan_full_story_2026-10-22.md` centers an MU/LUMIERE growth-forecast study and is not the current plan. The RANO-derived form and appendix remain unchanged.

## 5. Evidence and workflow rules

1. Cite every reported number to its output file or ledger entry.
2. Label results as observed, exploratory, simulated, assumed, or post hoc.
3. Keep negative results and state the weakest point beside each result.
4. Do not reuse the H-2 or H-3 evaluation sets as fresh confirmatory tests; those sets are spent.
5. Do not claim patient-specific treatment benefit, clinical validity, or a validated digital twin from the current simulations.
6. Before each new analysis, close the applicable source, provenance, inclusion, model, evaluation, and failure-rule gates in the current plan and its protocol manifest.
7. Keep patient data de-identified and follow each dataset's access and use terms.

## 6. Separate preserved work

The probabilistic RANO project is no longer the main project. Its plans, declarations, results, technical note, code, and outputs remain in place as a separate historical record. In particular, analysis_plan_v1.md and the RANO-derived research-plan drafts remain unchanged and do not define the new GBM project.

The old GBM source code, data, outputs, and archived documents remain in their existing locations. The current plan documents intended future analyses; it does not change historical results or authorize use of files whose terms remain unresolved.

## 7. Coding status and next work

- **Started:** Track A donor-split runner now has import-safe behavior and validates saved arrays, indices, labels, IDs, and patient disjointness before creating output. Synthetic tests pass. No real data were scored.
- **Shared split/provenance (2026-10-22):** `src/leakage_checks.py` now rejects blank, non-integer, negative, out-of-range and empty folds, splits that differ from the script-94 hash rule, and patients with no fold. Track B `src/100_pde_manifest.py` validates `split_mu.csv` on every load and reuses a cached forecast only when a SHA-256 signature of pair, fold, interval, treatment, grid or selected cell, solver and atlas matches; older unsigned caches are recomputed. Synthetic tests in `tests/test_split_provenance.py`; scratch-clone suite 142 passed. Script 100 was not run on patient data.
- **Track C simulator tests (2026-10-22):** `tests/test_twopop_simulator.py` checks the script-76 two-population model on synthetic parameters: exact reduction to the one-population model at zero resistance, irrelevance of the split with no cost and no drug or with eps = 1, drug immunity at full resistance, selection for resistance under drug, density bounds, first-order Euler convergence to the logistic solution (default 0.2-day step error about 2.5% at rho = 0.05 over 100 days), budget and window caps, paced spending and adaptive switching. Script 76 no longer creates its output folder on import. Scratch-clone suite 156 passed. These test the assumed model only, not drug effects.
- **Track A donor CV helper (2026-10-22):** new `src/donor_cv.py` resolves sample-to-donor identity with aliases and refuses conflicts; builds checked donor-grouped folds; runs nested donor-grouped selection (inner folds use outer-training donors only, one outer refit); and scores per donor. Synthetic tests in `tests/test_donor_cv.py` show no fit ever sees outer-test rows and a donor-lookup model falls to chance. Scratch-clone suite 162 passed. Not yet wired into any Track A script; no real data scored.
- **Track B synthetic parameter recovery (2026-10-22):** new `src/fk_mask_inverse.py` fits isotropic Fisher-Kolmogorov (rho, D) to a series of binary masks with the TensorFK solver. `tests/test_fk_parameter_recovery.py` makes reference masks at 1 mm, block-averages them to 2 mm, and fits at 2 mm (ellipsoid octant, scans at days 0/60/120/180). With three follow-ups the Jacobian has full rank, and equal-front-speed pairs give different scans, so rho and D were separately identifiable in this setup. With one follow-up a ridge remains: only a combination is constrained. Fits are off the bounds, at least as good as the truth, Dice >= 0.95. Log errors (rho, D): 0.03, 0.066 at truth (0.03, 0.1); 0.162, 0.408 at truth (0.05, 0.03), where D is about 33% low. Both stay within the linearised bias from the measured 1 mm vs 2 mm mismatch (up to 5% log volume at the truth). Weakest point: 2 mm discretisation alone biases D by up to about 1.5x in the second case. The script-51 volume surrogate puts D on its lower bound and misses rho by more than 2x on the same PDE volumes (kept as a limitation test). Scratch-clone suite 172 passed (the recovery test takes about 2 minutes). No patient scans fitted; this is not evidence of recovery from real scans.
- **Next Track B step:** run a prespecified synthetic scenario grid with clean/noisy training masks and isotropic/anisotropic generating truth. Compare PDE, persistence, and volume-trend forecasts on the same later held-out mask. Compare calibrated error bounds from ensemble boundary disagreement against training spatial and volume residual scores using separate calibration and evaluation cases. Report coverage, bound width, forecast losses, and failure cases. Do not fit patient scans until this synthetic check and the image-level dataset manifest gate pass.
- **Track B synthetic forecast and uncertainty helpers (2026-10-23):** `src/fk_mask_inverse.py` has a tensor-field forward helper for anisotropic synthetic references. `src/fk_synthetic_validation.py` adds reproducible signed-distance boundary perturbations, a deterministic volume-trend baseline, held-out next-scan comparison, bootstrap forecast ensembles, a spatial training-residual score, and split-conformal upper bounds for normalized forecast error. The ensemble and calibration helpers use training/calibration cases only; held-out targets enter only forecast scoring. Boundary noise is an assumed simulation model, not measured segmentation error. Scratch-copy checks passed: 24 focused tests, including parameter recovery. No patient data were read or fitted, and no synthetic scenario-grid result or calibrated-coverage result has been measured yet.
- **Spent evaluation:** Do not rerun `src/20_patient_split_methods.py` on the saved five-patient test split and call it fresh confirmation. It is historical/exploratory only.
- **Next shared step:** inspect and stabilize data manifests, source-file identities, split hashes, run configuration, and output provenance before new model runs.
- **Track A:** resolve donor/sample aliases and dataset terms; implement nested donor-grouped selection, donor-level evaluation, panel mapping, and held-out spatial transfer without label or identity leakage.
- **Track B:** reconcile patient/visit manifests; add solver and physical-unit tests; compare PDE, volume, and recurrence models with fixed baselines; test parameter recovery and forecastability before patient fitting.
- **Track C:** freeze simulator equations, seeds, budgets, and endpoints; validate two-population dynamics; compare adaptive, fixed, and learned policies on common simulated cases; report sensitivity and failures. Keep MGMT association observational and separate.
- **Validation:** use synthetic/unit tests first. Run full experiments only after input rights, provenance, inclusion, split, and failure rules are recorded. Preserve old outputs; never reuse spent test sets as confirmatory.

The research-plan form and appendix specify intended methods, data gates, primary measures, baselines, and failure rules. Keep later code aligned with them; revise the plan before changing the research design. Do not claim clinical benefit, validated patient-specific predictions, or a fused digital twin.
