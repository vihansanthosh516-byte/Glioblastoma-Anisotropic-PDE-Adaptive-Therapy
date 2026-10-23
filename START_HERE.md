# START HERE

Updated 2026-10-22. Read this file before continuing work. The main project is now the earlier GBM Tracks A, B, and C work for a local science fair.

## 1. Current goal

Reorganize and improve the old GBM project into a defensible science-fair research plan. The future plan must account for Tracks A, B, and C, use openly downloadable data, and clearly separate measured data from assumptions and simulations.

This update only organizes files and records known improvement needs. It does not choose a final research question, approve an analysis, run code, download data, or change results.

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

Prior drafts exist. The forward plan in paper/research_plan_full_story_2026-10-22.md centers a MU/LUMIERE growth-forecast study and does not set roles for all three tracks, so it is not adopted as the current plan. The RANO-derived form and appendix remain unchanged. Use the improvement backlog when preparing the new GBM plan.

## 5. Evidence and workflow rules

1. Cite every reported number to its output file or ledger entry.
2. Label results as observed, exploratory, simulated, assumed, or post hoc.
3. Keep negative results and state the weakest point beside each result.
4. Do not reuse the H-2 or H-3 evaluation sets as fresh confirmatory tests; those sets are spent.
5. Do not claim patient-specific treatment benefit, clinical validity, or a validated digital twin from the current simulations.
6. Before any new analysis, write a separate GBM plan that fixes the question, open dataset and license, inclusion rules, model, baselines, held-out evaluation, and failure rule.
7. Keep patient data de-identified and follow each dataset's access and use terms.

## 6. Separate preserved work

The probabilistic RANO project is no longer the main project. Its plans, declarations, results, technical note, code, and outputs remain in place as a separate historical record. In particular, analysis_plan_v1.md and the RANO-derived research-plan drafts remain unchanged and do not define the new GBM project.

The old GBM source code, data, outputs, and archived documents also remain in their existing locations. This organization pass changes navigation and project status only.

## 7. Next step

Use the improvement backlog to build the GBM science-fair research plan. Decide which track supplies the primary question and outcome, and define the roles of the other tracks. Do not start new analysis or code until that plan is reviewed and recorded.
