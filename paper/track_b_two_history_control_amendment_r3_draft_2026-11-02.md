# Track B two-history control amendment, r3 draft (PROPOSED; 2026-11-02)

**Status: PROPOSED, not adopted.** This revision addresses the remaining score-definition and sample-size gaps. It authorizes no patient tests, jobs, patient-data reads, downloads, contacts, or scoring. It is not a cohort-specific protocol. It supersedes r2, which did not name the candidate score used by the max-T family or set a minimum inner-training size.

The earlier r1 draft [`track_b_two_history_control_amendment_draft_2026-11-02.md`](track_b_two_history_control_amendment_draft_2026-11-02.md) is also superseded: its cross-outer-fold fallback can leak held-out outcomes, its volume fallback is undefined, and it mixes Duke and LUMIERE rules. The synthetic control remains unadopted and is not wired to a patient runner. Any future evaluation requires a dated amendment for one selected cohort, including its label map, support, eligibility rules, folds, comparator family, source gates, compute gate, and authorization.

## 1. Question and claim

The proposed primary comparison is between:

- **Candidate:** PDE-derived spatial ranking.
- **Control:** a two-history signed-distance shape extrapolation.

Both use the same history-only volume forecast draws, origin support, patient set, target labels, and forecast interval. The focused question is whether the candidate’s ranking gives higher next-observed-visit core Dice than the specified control at those shared volume draws.

This compares predictive systems. It does **not** identify a causal PDE effect: the systems differ in fitted structure, and any boundary-perturbed analysis uses model-specific uncertainty. The result is agreement with the cohort’s source labels, not true tumor, clinical benefit, generalization, or novelty.

If the realized target-visit interval is supplied, the estimand is retrospective and conditional on that interval. It is not a prospective fixed-day forecast.

## 2. Common volume draws and failure policy

- Fit the volume model using training patients only. Within inner selection, refit it using inner-training patients only.
- Inputs are histories and their intervals plus the realized forecast interval. No target mask, target volume, or held-out outcome enters a forecast or fit.
- Generate a patient’s volume draws once, then pass the identical draws to every arm that requires them. Match each draw inside common support S. The protocol must record the seed and the requested voxel counts.
- A history-only failure of the shared volume forecast (nonfinite, nonpositive, or outside the valid support capacity) invokes one common, predeclared fallback for all volume-matched mask arms: the visit-2 core clipped to S. Keep the patient in the denominator and record the stage and exact reason.
- An empty or full visit-1 core is detected before target access. For the primary pair, use the same visit-2 persistence fallback for both arms, keep the patient in the denominator, and count the event. Do not let the control use f = 0 while the candidate separately fails.
- For the focused candidate-versus-control pair, a prespecified recoverable failure in either shape ranking also gives both primary arms that same persistence mask for that patient. This is a conservative common-failover system comparison; report its frequency and reasons. Do not silently drop the patient or call that row a successful shape forecast.
- Catch only named, expected numerical/model failures. Unexpected exceptions, data-integrity failures, or code errors stop the run as incomplete; they do not trigger fallback.
- Empty or full visit-2 support is an eligibility issue, not a runtime tuning choice. Decide the rule before target access, count affected patients, and apply the minimum before any target mask is opened.

## 3. Nested selection without leakage

- The proposed control grid is `SHRINK_TAU_DAYS {0, 7, 28} × MAX_FACTOR {1, 4, 8} × MAX_SD_DELTA_MM {10, 30}`, plus f = 0, the last-mask signed-distance rule.
- This grid was designed after the synthetic sensitivity screen was inspected. The screen’s winning cell is not carried forward. The grid and selection rule must be frozen before any patient outcome is scored; do not describe them as chosen before all data or before synthetic development.
- For each outer fold, select settings using only that fold’s outer-training patients, with patient-grouped inner folds. Refit the volume model separately inside each inner training split. Score each setting by inner-validation patient-macro pre-boundary Dice.
- Select the setting with the highest mean inner-validation Dice. If values differ by at most 0.005, choose the more damped setting using this order: f = 0, lower cap, higher tau, lower clip.
- The outer-test patient is scored once with the setting selected without that patient.
- **Minimum inner-training size (proposed):** every inner-training split must contain at least 10 eligible patients. If any split is smaller, stop the full analysis before opening target masks. This is a feasibility floor carried from the earlier proposal, not a power or precision guarantee. No leave-one-out exception is allowed.
- Across-fold setting frequencies and score margins are descriptive only. **No pooled outer-fold stability statistic may alter any outer-fold prediction.** In particular, remove the prior rule that fell back across every fold when fewer than three outer folds picked the same setting; it leaks because a test patient is training data for other outer folds.
- If the predeclared minimum inner-training size cannot be met in any fold, stop the full analysis before target scoring. Do not create fold-specific exceptions after outcomes are seen.

## 4. Scores, baselines, and decision rule

- **Focused primary contrast:** candidate minus two-history control, patient-level macro Dice, pre-boundary.
- **Exact candidate score for the primary family:** threshold the candidate's `pre_boundary_probability` at inclusion probability >= 0.5 inside S. The resulting binary mask is `candidate_pre_boundary`. Use this same candidate score for the focused contrast and every adjusted max-T contrast. The runner's current `mask_draws` score includes boundary perturbation; it is a different candidate score and must not enter this family.
- **Separate boundary score:** threshold the mean of the post-boundary candidate `mask_draws` at 0.5 inside S. Report it as a secondary result only. Never mix a pre-boundary candidate score and post-boundary candidate score in one max-T family.
- Report the boundary-perturbed version separately as secondary; do not pool it with pre-boundary Dice. Boundary results are not a neutral isolation of shape because covariance may be PDE-calibrated.
- Report oracle-volume shape diagnostics in a separate descriptive table. They are not forecasts and cannot support the primary claim.
- The cohort-specific amendment must name every required mask comparator and its exact recipe. Do not import LUMIERE’s `COMPARATORS` constant or Duke label rules into another cohort by reference.
- **Proposed strict pass gate:** preserve the existing project rule that the candidate must beat every required mask baseline. Include the focused control and every other required comparator in one predeclared max-T family, all compared with `candidate_pre_boundary`. The run passes only if every adjusted lower confidence bound is above zero; any missing arm or nonfinite result makes the analysis incomplete. If the authors instead want a primary-only success rule, that is a separate amendment and must narrow the claim to the primary contrast.
- Count prediction remains on its count metrics only; it is not converted to a Dice mask.

## 5. Cohort separation

This document defines only a candidate control concept. It does not define a runnable Duke or LUMIERE protocol.

- **Duke:** no written access grant is recorded. Any future Duke protocol must resolve access and derived-output rights, school approval, prior-use overlap, eligible three-visit patient count, treatment/interval policy, label meanings, support, registration QC, folds, and run time. Do not read restricted data until its own gates and user authorization pass.
- **LUMIERE:** it is previously exposed and can support exploratory re-analysis only. Its core labels and source uncertainties differ from Duke’s. Do not call it independent confirmation. Do not reuse Duke’s 1|3 core/label-4 cavity rule; LUMIERE’s recorded target rule is labels 1|2 with no documented cavity label.
- **History eligibility:** empty/full visit 1 uses the common primary-pair fallback above. Empty/full visit 2 within S is excluded before target access because it cannot define the proposed shape ranking. Report each count and stop if the remaining cohort fails its prespecified minimum.
- The chosen cohort’s source labels are the endpoint. The protocol must state what label uncertainty means for interpretation.

## 6. Code and validation required before adoption

1. **Settings API:** implemented in the synthetic prototype as `TwoHistoryShapeSettings` with finite/range validation and JSON-friendly recording. Former defaults are preserved. The development screen was synthetic; it does not validate patient performance or select the winning cell.
2. **Fold-local selector: implemented as a pure synthetic utility.** `validate_inner_fold_plan` checks exact training/validation complements, one validation appearance per outer-training patient, no outside IDs, and the 10-patient training floor. `select_two_history_settings` accepts only validation scores for IDs in that fold, aggregates patient-macro Dice, and applies the 0.005 tie rule. Outer-test IDs are rejected. It does not load targets, fit models, or wire into a patient runner.
3. **Shared history-only volume forecast: implemented.** `history_only_volume_forecast` fits and samples volume draws without a PDE call. The candidate accepts the precomputed object and validates its history, interval, fit, seed, draw count, and support capacity signature. An output fingerprint detects changed volume counts or derived values. Synthetic tests confirm the existing candidate seed path gives identical requested volumes and forecast outputs. Common paired-arm fallback behavior is not wired.
4. **Typed expected errors: implemented.** Capacity, invalid-volume, and unusable-boundary failures have named exception types. The LUMIERE PDE provider converts only these expected failures to `ForecastFailure`; unexpected `ValueError` propagates. This does not yet implement the proposed common fallback for both shape arms.
5. Add pre-target eligibility for the selected cohort, including empty/full visit-1 and visit-2 handling and the minimum-patient stop.
6. Add `candidate_pre_boundary` to the future result contract and use it as the candidate key for the primary max-T family. Store the post-boundary candidate score separately as secondary. Do not wire this proposal into a patient runner until a cohort amendment is adopted.
7. Reconcile the result writer with the exact all-baseline pass rule. Record per-patient fallbacks, missing arms, the candidate score version, and selected control cell.
8. Before adoption, calculate expected interval precision at the expected eligible cohort size using only development data and patient-level paired differences. State the chosen precision target and assumptions. If no defensible precision target or eligible count exists, label the work feasibility/exploratory; do not claim it is powered.
9. Recompute run time with representative history and forecast gaps, all inner calibration calls, all required models, refits, prior settings, and checkpoint/restart overhead. Use wall-clock hours on the named PC; state the scenario assumptions.

**Synthetic verification update (2026-11-02):** five focused suites passed in a scratch copy: 68 passed across the hybrid forecast, LUMIERE-provider, two-history, MU pipeline-component, and patient-CV suites. These tests use synthetic fixtures. No patient runner was changed, no patient data was read, and no score was produced.

No patient test is authorized until source, access, school, cohort, support, registration, protocol-freeze, input-hash, and user-authorization gates pass for one named cohort.

## 7. Precision and power gate

The Duke eligible-patient count is unknown. Before adoption, report expected confidence-interval width for the paired patient-level contrasts across plausible eligible sample sizes, using only prior development summaries or prespecified variance scenarios. Do not use future test outcomes. Set a precision target before access to target labels. A minimum training size of 10 does not establish adequate power. If the expected precision is inadequate, stop the confirmatory framing or label the work a feasibility study.

## 8. Required interpretation limits

Even if the proposed comparison passes, the report may say only that the candidate system outperformed the named baselines on the selected cohort’s source labels and observed interval distribution. It may not say PDE caused the gain, that the labels are true tumor, that the method is clinically useful, that it generalizes, or that the approach is novel without a separate current prior-art analysis.
