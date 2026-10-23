# GBM Tracks A/B/C improvement backlog

**Review date:** 2026-10-22

**Scope:** Repository records and audits only. No web research, code, data, or model runs were used for this note. This is input for a future science-fair plan, not a new protocol or novelty claim.

## Track A — Single-cell and spatial data

**Documented weakness.** The headline cell-level result came from a random-cell split that leaked patients: reported accuracy was 78.7%, while a patient-ID-only model reached 82.8%. The patient-level C-GAT result was 63.8%, close to the cells-per-patient baseline at 64.6%; the audit flags donor composition as a confound. The reported atlas has 25 patients, and no matched patient link from this atlas to the MRI cohorts is established. Age/survival results are observational associations; the audit says age is already known and the adjusted external result is not predictive validation. MGMT analysis was underpowered.

**Evidence:** `README_STATUS.md` sections 2–3; `paper/codex_audit/old_tracks.md` findings on C-GAT and age; `vault/PAPER_DATA_SOURCES.md`; `archive_old_gbm/TRACK_A_RESULTS.md`.

**Upgrade ideas for a new plan.** Make the patient, not the cell, the unit of generalization. Use patient-grouped folds, compare against patient-ID and cell-composition baselines, balance or reweight donors, and report variation across held-out patients. Define one target that the atlas actually labels. Keep expression-based survival association separate from prediction. Treat any molecular-to-MRI connection as a population-level hypothesis unless shared patient identifiers and provenance are proved.

**Still blocked.** Only 25 atlas patients are listed in the data-source record. The old notes do not establish an external cohort with the same spatial labels, a reliable patient-level clinical endpoint for the atlas, or a valid cross-track patient join. Exact dataset release terms and local file provenance need checking before reuse.

## Track B — Tumor growth and imaging

**Documented weakness.** Historical PDE growth forecasts did not show an advantage over simple baselines. The GBM-only PDE-versus-persistence estimate was +0.013 Dice with an interval crossing zero; the geometric-baseline comparison was also inconclusive. LUMIERE growth-only forecasts were worse than no-change in the recorded analysis. In the fitted-volume inverse problem, diffusion hit a bound in 147/154 cases and growth rate hit its lower bound in 124/154; the notes say diffusion was not identifiable in that setup. Some cohort counts differ by workflow, and the exact inclusion lineage is incomplete.

**Evidence:** `README_STATUS.md` sections 2–3; `vault/PAPER_DATA_SOURCES.md`; `archive_old_gbm/vault/TRACK_B_RESULTS.md`; `archive_old_gbm/paper/summary_all_tracks.md`; `vault/PAPER_FINDINGS_LEDGER.md` entries cited by those summaries.

**Upgrade ideas for a new plan.** Start with an auditable patient-by-visit table: dates, masks, sequences, and outcome must align. Require enough visits for the declared forecast; report exclusions. Compare every new model with persistence and a simple geometric or population-rate baseline on the same held-out patients and forecast times. Test parameter recovery on simulated data before fitting patient data. If parameters remain weakly identified, target uncertainty in the future forecast itself rather than claiming patient-specific biological parameter estimates. The archived trusted-twin draft proposes simulation checks, forecastability/abstention ideas, and explicit failure reporting; these are proposals only.

**Still blocked.** Old test sets H-2/H-3 are spent and cannot serve as fresh locked tests. The repository does not establish a new independent test cohort with a verified patient-level count of sufficiently complete visits. A new dataset’s license, keys, segmentation quality, and eligibility must be checked before choosing it. Numerical solver checks do not validate biological predictions.

## Track C — Treatment suppression and adaptive control

**Documented weakness.** Treatment effects and kill rates are assumed in the simulations. The `real_test` policy run uses fitted growth rates for 21 patients but the same assumed kill-rate table; it is not a real-drug evaluation. The policy did not win on time-to-progression, and a simple rule matched its day-90 volume result. The adaptive-therapy simulation reduced modeled resistance in some comparisons while producing higher final tumor mass and earlier progression for some patients. The model lacks clonal competition, and Track C parameter fits inherit Track B identifiability limits.

**Evidence:** `README_STATUS.md` sections 2–3; `archive_old_gbm/vault/TRACK_B_RESULTS.md`; `archive_old_gbm/vault/TRACK_C_RESULTS.md`; `paper/codex_audit/old_tracks.md`; `vault/claim_ledger.md`.

**Upgrade ideas for a new plan.** Keep treatment effects explicitly simulated unless patient-level dose, dates, and response are available. Compare schedules under the same total dose and stated toxicity assumptions. Include fixed-treatment and simple-rule baselines. Use time-to-progression and tumor burden together, and report cases where the policy harms either outcome. Add uncertainty or sensitivity over kill rates before interpreting policy rankings. A two-population sensitive/resistant model is an existing repository direction, but it still needs empirical parameter support or must remain a simulation.

**Still blocked.** No cited repository source establishes a fully open cohort with joined patient-level MRI, treatment dose/timing, and response adequate to estimate individual drug effects. No result supports a clinical recommendation, patient-specific treatment benefit, or biological biomarker cutoff.

## Cross-track project design

The records do not support one fused patient-level model. Track A’s atlas patients are not linked to the MRI patients; the recorded A-to-B score is constant at 1.0 for all 61 cases, and its source lineage is unresolved (`vault/PAPER_DATA_SOURCES.md`; `paper/codex_audit/old_tracks.md`). A new plan can preserve A, B, and C as related GBM workstreams, but must choose one primary question and state which track supplies its primary data and endpoint. Other tracks can be labeled exploratory or simulation-only. Do not call the combined system a validated digital twin.

The archived `project_plan_trusted_twin.md` is a historical proposal. Its theorem and novelty ideas were not established. `README_STATUS.md` and `paper/codex_audit/old_tracks.md` are the better starting point for historical outcomes and corrections. Old documents with “first,” “novel,” “revolutionary,” or “validated” claims should remain records, not active project claims.

Before any new analysis, the science-fair plan still needs: one primary research question; an exact public dataset and license; patient-level inclusion rules; a fixed mathematical model and baselines; a held-out evaluation that avoids spent test sets; a declared failure criterion; and local fair review of secondary-data and AI-use requirements. No new code or analysis is approved by this backlog.
