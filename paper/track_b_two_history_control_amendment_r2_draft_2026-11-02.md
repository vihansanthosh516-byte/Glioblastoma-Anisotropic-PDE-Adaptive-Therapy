# Track B two-history control amendment, r2 draft (PROPOSED; 2026-11-02)

**Status: PROPOSED, not adopted.** This is a technical correction to the first proposal. It authorizes no tests, jobs, patient-data reads, downloads, contacts, or scoring. It is not a cohort-specific protocol. The first draft [`track_b_two_history_control_amendment_draft_2026-11-02.md`](track_b_two_history_control_amendment_draft_2026-11-02.md) is superseded for future consideration because its cross-outer-fold stability fallback can leak held-out outcomes and its cohort-specific rules are mixed.

This revision is based on a code/science audit. The existing synthetic control remains unadopted and is not wired to a patient runner. Any future evaluation requires a dated amendment for the selected cohort, including its own label map, support, patient minimum, folds, comparator family, source gates, and authorization.

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
- For the focused candidate-versus-control pair, a prespecified recoverable failure in either shape ranking also gives both primary arms that same persistence mask for that patient. This is a conservative common-failover system comparison; report its frequency and reasons. Do not silently drop the patient or call that row a successful shape forecast.
- Catch only named, expected numerical/model failures. Unexpected exceptions, data-integrity failures, or code errors stop the run as incomplete; they do not trigger fallback.
- Empty or full visit-2 support is an eligibility issue, not a runtime tuning choice. Decide the rule before target access, count affected patients, and apply the minimum before any target mask is opened.

## 3. Nested selection without leakage

- The proposed control grid is `SHRINK_TAU_DAYS {0, 7, 28} × MAX_FACTOR {1, 4, 8} × MAX_SD_DELTA_MM {10, 30}`, plus f = 0, the last-mask signed-distance rule.
- This grid was designed after the synthetic sensitivity screen was inspected. The screen’s winning cell is not carried forward. The grid and selection rule must be frozen before any patient outcome is scored; do not describe them as chosen before all data or before synthetic development.
- For each outer fold, select settings using only that fold’s outer-training patients, with patient-grouped inner folds. Refit the volume model separately inside each inner training split. Score each setting by inner-validation patient-macro pre-boundary Dice.
- Select the setting with the highest mean inner-validation Dice. If values differ by at most 0.005, choose the more damped setting using this order: f = 0, lower cap, higher tau, lower clip.
- The outer-test patient is scored once with the setting selected without that patient.
- Across-fold setting frequencies and score margins are descriptive only. **No pooled outer-fold stability statistic may alter any outer-fold prediction.** In particular, remove the prior rule that fell back across every fold when fewer than three outer folds picked the same setting; it leaks because a test patient is training data for other outer folds.
- If the predeclared minimum inner-training size cannot be met in any fold, stop the full analysis before target scoring. Do not create fold-specific exceptions after outcomes are seen.

## 4. Scores, baselines, and decision rule

- **Focused primary contrast:** candidate minus two-history control, patient-level macro Dice, pre-boundary.
- Report the boundary-perturbed version separately as secondary; do not pool it with pre-boundary Dice. Boundary results are not a neutral isolation of shape because covariance may be PDE-calibrated.
- Report oracle-volume shape diagnostics in a separate descriptive table. They are not forecasts and cannot support the primary claim.
- The cohort-specific amendment must name every required mask comparator and its exact recipe. Do not import LUMIERE’s `COMPARATORS` constant or Duke label rules into another cohort by reference.
- **Proposed strict pass gate:** preserve the existing project rule that the candidate must beat every required mask baseline. Include the focused control and every other required comparator in one predeclared max-T family. The run passes only if every adjusted lower confidence bound is above zero; any missing arm or nonfinite result makes the analysis incomplete. If the authors instead want a primary-only success rule, that is a separate amendment and must narrow the claim to the primary contrast.
- Count prediction remains on its count metrics only; it is not converted to a Dice mask.

## 5. Cohort separation

This document defines only a candidate control concept. It does not define a runnable Duke or LUMIERE protocol.

- **Duke:** no written access grant is recorded. Any future Duke protocol must resolve access and derived-output rights, school approval, prior-use overlap, eligible three-visit patient count, treatment/interval policy, label meanings, support, registration QC, folds, and run time. Do not read restricted data until its own gates and user authorization pass.
- **LUMIERE:** it is previously exposed and can support exploratory re-analysis only. Its core labels and source uncertainties differ from Duke’s. Do not call it independent confirmation. Do not reuse Duke’s 1|3 core/label-4 cavity rule; LUMIERE’s recorded target rule is labels 1|2 with no documented cavity label.
- The chosen cohort’s source labels are the endpoint. The protocol must state what label uncertainty means for interpretation.

## 6. Code and validation required before adoption

1. **Settings API:** implemented in the synthetic prototype as `TwoHistoryShapeSettings` with finite/range validation and JSON-friendly recording. The former defaults are preserved. Twenty-one focused synthetic tests pass in a scratch copy. This does not select settings on patient data or wire the control into a runner.
2. Implement a fold-local selector that accepts only outer-training IDs and returns the selected cell plus inner results. Add a leakage test that changes outer-test targets and proves selection and forecasts do not change.
3. Factor the shared history-only volume forecast so the candidate and control consume the exact same draws, including fallback cases.
4. Implement the common failure rules above with typed expected failures. Test expected fallback and unexpected abort behavior separately.
5. Add pre-target eligibility for the selected cohort, including empty/full origin-core handling and the minimum-patient stop.
6. Wire the control only into the selected cohort runner. Test identical patients, folds, support, target indexing, volume draws, and comparator completeness on synthetic data in a scratch copy.
7. Reconcile the max-T family and result writer with the exact all-baseline pass rule. Record per-patient fallbacks, missing arms, and the selected cell.
8. Recompute run time with representative history and forecast gaps, all inner calibration calls, all required models, refits, prior settings, and checkpoint/restart overhead. Use wall-clock hours on the named PC; state the scenario assumptions.

No patient test is authorized until source, access, school, cohort, support, registration, protocol-freeze, input-hash, and user-authorization gates pass for one named cohort.

## 7. Required interpretation limits

Even if the proposed comparison passes, the report may say only that the candidate system outperformed the named baselines on the selected cohort’s source labels and observed interval distribution. It may not say PDE caused the gain, that the labels are true tumor, that the method is clinically useful, that it generalizes, or that the approach is novel without a separate current prior-art analysis.
