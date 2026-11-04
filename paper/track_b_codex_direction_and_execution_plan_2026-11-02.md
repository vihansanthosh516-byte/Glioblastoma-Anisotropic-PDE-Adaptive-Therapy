# Track B direction and execution plan (Codex review, 2026-11-02)

**Status: decision plan only.** This document does not adopt an amendment, authorize patient-data reads, authorize a run, or change a result. Track A and Track C are outside scope. No patient voxels were read for this review.

## Track B goal

Test whether a reaction-diffusion/PDE model contributes useful **spatial tumor-shape forecasting** beyond strong history-based methods, when every method gets the same patient history, support, forecast interval, and training-only volume forecast. The result must be reproducible, patient-level, and reported even when it is negative.

The claim is about prediction of source-provided labels. It is not a causal claim about PDE biology, a claim about true tumor, a clinical claim, or a broad novelty claim.

## Current scientific position

- Broad “PDE plus learned correction” novelty is unavailable: [SoC-DT](https://arxiv.org/html/2510.03287) already evaluates a PDE/neural-correction approach on real patients. Track B must make a narrower, measured contribution or present itself as a replication/engineering study.
- The latest proposed narrow question is useful: at identical history-only volume draws, does PDE-based spatial ranking outperform a two-history non-PDE shape extrapolation? This tests two predictive systems; it does not isolate a causal PDE effect.
- The primary target should be **next observed visit conditional on its realized interval** unless the cohort and visit schedule support a fixed, origin-known horizon. Call it retrospective horizon-conditioned evaluation, not a prospective fixed-day forecast.
- The existing LUMIERE evidence is exploratory and exposed. The project previously used LUMIERE, and LUMIERE appears in [PREDICT-GBM](https://arxiv.org/html/2509.13360). Its automated labels, unresolved label-0 meaning, undocumented tool version, unresolved support/registration rules, and long runtime weaken the value of a large new patient run. Its public [Figshare record](https://springernature.figshare.com/articles/dataset/LUMIERE_dataset_-_MRI_data_and_automated_segmentations/21249516) also needs reconciliation with the project’s existing non-commercial/derived-output restrictions before any reuse or sharing.
- Duke is the strongest structural match found so far: the [official record](https://zenodo.org/records/20548528) lists 61 patients and 3–18 sessions with masks. **Status update (2026-11-03):** the user reports that download access was granted and download started. At the latest local check, `C:\Users\vihan\DATASET_FINAL.zip` was 66,347,008 bytes, below the expected 17,043,569,844-byte archive size, and the file was growing. The public record still marks files restricted and describes two license paths, CC BY-NC-ND 4.0 or a custom Duke license; the accepted terms are not yet saved. The authors' [usage guide](https://github.com/deshanyang/GBM-DIR-QA) says visits are not registered to each other and relative dates may be offset by days. Triplet count, patient overlap, derived-output terms, school approval, and longitudinal MRI registration QC remain unresolved. The user authorized metadata-only local inventory after download completion; MRI, mask, and landmark-coordinate contents remain outside this scope. See [Duke readiness note](track_b_duke_readiness_2026-11-03.md). Duke remains **not ready for patient scoring**.
- A separate public-record screen found no verified fresh, no-request cohort that passes the combined requirements for at least three masked visits, clear reuse terms, GBM eligibility, and project independence. [RHUH-GBM](https://www.cancerimagingarchive.net/collection/rhuh-gbm/) is structurally close but is already exposed through PREDICT-GBM. BraTS Longitudinal and WAW-Glioma remain unverified on access/terms, eligible GBM counts, longitudinal masks, or overlap. This was a bounded search, not proof that no such dataset exists.

## Recommended direction

**Do not run the full four-analysis LUMIERE design.** It would add an expensive, previously exposed, automated-label result while source and protocol gates remain open. The prior runtime figures are projections, not guarantees, and the long-gap PDE cases can dominate them.

Keep one conditional external path: **Duke only if written access and terms are granted and every pre-run gate passes.** First resolve eligibility and project overlap from permitted metadata. If Duke fails a gate, close the external-patient branch and report external confirmation as unavailable. Do not continue open-ended cohort scouting or switch to an already exposed cohort to obtain a positive result.

If no independent cohort becomes available, the defensible endpoint is a methods/simulation report plus transparent prior exploratory results. Do not describe synthetic success as patient validation or as high novelty. More synthetic tuning cannot repair the absence of independent patient evidence.

## Corrections required before the two-history proposal can be implemented

The first two proposed amendments are **not safe to adopt as written**. The current proposal is [`track_b_two_history_control_amendment_r3_draft_2026-11-02.md`](track_b_two_history_control_amendment_r3_draft_2026-11-02.md). It addresses:

1. **Outer-fold leakage:** a global rule based on how many outer folds choose a setting lets each test patient affect other folds’ selections and then change that patient’s own prediction. Stability can be descriptive only, or computed entirely inside each outer-training set. It must never trigger a cross-fold fallback.
2. **Undefined failures:** name the exact history-only fallback, the expected recoverable errors, what is counted, and which errors abort the run. Do not silently catch programming errors. Keep every eligible patient in the denominator and report fallbacks.
3. **Cohort mismatch:** the draft mixes Duke labels (core 1|3, cavity 4) with the LUMIERE runner and its comparator family (LUMIERE core 1|2; no cavity label). The control proposal is now cohort-neutral. Any run needs a separate cohort-specific amendment, label map, runner, and analysis family.
4. **One coherent decision rule:** define whether success requires the focused candidate-vs-control contrast only, or all required mask baselines. The proposed revision keeps the project’s strict all-required-baselines gate and names the complete max-T family; if the family is incomplete, the result is incomplete.
5. **Development disclosure:** the settings grid was designed after a synthetic screen had been inspected. The revised text says this plainly. It does not carry forward the screen’s winning setting.
6. **Eligibility before target access:** require valid, nonempty, non-full origin support for the shape comparison, or define a justified alternative before adoption. Count exclusions and apply the minimum before opening targets.
7. **Claim boundary:** the endpoint is agreement with automated/source labels at the next observed visit conditional on interval, not truth, causal PDE effect, clinical value, or prospective fixed-day performance.
8. **One candidate score:** the proposed max-T family uses the same pre-boundary candidate mask (`pre_boundary_probability >= 0.5` inside S) for every contrast; the current post-boundary candidate score is separate and secondary.
9. **Feasibility:** a proposed minimum inner-training size of 10 is explicit, empty/full history handling is stated, and a patient-level precision check is required before adoption. The 10-patient floor is not a power claim.

## Ordered work and stop/go gates

### Stage 1 — Freeze the scientific question on paper

- Resolve the corrected two-history proposal and its relation to the existing all-baseline pass rule.
- Pick one cohort and write its own protocol. Do not blend Duke and LUMIERE labels, support, time windows, comparators, sample minimums, or compute assumptions.
- Freeze the target definition, primary metric, patient-level aggregation, uncertainty method, folds, missingness/fallback rules, multiplicity family, and what counts as a pass.
- Keep all existing patient results and the synthetic screen visible in the prior-result/exposure ledger. No screen result is a hidden confirmatory result.

**Stop here** if the scientific question cannot be stated without a causal or novelty claim the design cannot test.

### Stage 2 — Rights, approval, and cohort eligibility

For Duke, before any restricted-file access:

- Save written access approval and the exact license/agreement, including allowed derived masks, predictions, model weights, figures, and summary outputs.
- Save written school approval with date, approver, and conditions.
- Resolve patient-level overlap with project-used cohorts, within the permitted access scope.
- Establish the number of independent patients with three usable visits, eligible labels, usable intervals, and acceptable metadata. Set a minimum before looking at target outcomes; justify it with precision/power or label it as a feasibility threshold.
- Resolve source-backed label meanings, valid support, treatment/surgery policy, and registration evidence.

If written access, output rights, approval, or the minimum eligible cohort fails, **stop the patient branch**. Do not substitute LUMIERE and call it independent.

### Stage 3 — Cohort-specific technical preflight

Only after Stage 2 and a separate user-approved read scope:

- Freeze image-only registration/support QC before target-label reads.
- Validate file geometry and selected inputs. Record exclusion counts, reasons, and hashes only for files the user has authorized to read.
- Freeze the eligible patient list and its hash before target scoring.
- If target support, image registration, or label definitions fail, stop; do not tune thresholds on target masks.

### Stage 4 — Synthetic implementation and compute gate

- **Done, synthetic utilities only:** `src/fk_two_history_shape.py` now has a validated settings object, an outer-training-only fold validator, and a nested selector that rejects IDs outside each validation fold. `src/fk_hybrid_forecast.py` now factors the history-only volume forecast and accepts its checked draws. Forecast-volume and unusable-boundary failures have typed errors; the LUMIERE provider converts only named expected failures.
- **Scratch verification:** five focused hybrid forecast, LUMIERE-provider, two-history, MU pipeline-component, and patient-CV suites passed: 68 tests. They use synthetic fixtures. The selector and shared-volume helper are not wired into a patient runner, and no patient data or scores were used.
- **Still required before any patient-run integration:** define common paired-arm fallbacks and pre-target eligibility in a cohort-specific protocol; add `candidate_pre_boundary` to a future result contract; reconcile the result writer; and test the full synthetic candidate-control path and outer-fold isolation.
- Run no patient-derived timing job until source and access gates pass. Use synthetic shapes/gaps representative of the permitted cohort, including the longest interval, and count every PDE calibration, forecast, inner fit, U-Net training/refit, and prior run.
- Set a practical wall-clock ceiling before the timing run. If the estimate exceeds it, stop and propose a scope amendment; do not weaken the strongest baseline after seeing outcomes.

### Stage 5 — Seal and one locked evaluation

- Freeze code, settings, folds, cohort list, package/runtime record, and input hashes at a clean commit before scoring.
- Run the frozen cohort once. Do not tune on outer-test outcomes or repeatedly inspect test results.
- Report every required comparator, donor/patient-level paired differences and intervals, all exclusions/fallbacks, model failures, runtime, and negative findings. Preserve results unchanged.
- If the full required arm family is not available, label the analysis incomplete. Never turn an incomplete result into a pass.

## Practical success definition

Track B can count as a successful **bounded research study** if one lawful, eligible cohort supports one locked patient-level evaluation; the candidate beats every required baseline under the frozen rule; code and inputs are sealed; and the conclusion is limited to that cohort’s source labels and visit intervals. A null or negative result can still be a valid completed study.

It does not become independently confirmatory, clinically useful, causal, or highly novel just because it runs. Those are separate claims requiring separate evidence.

## Current state

- Current LUMIERE and Duke patient scoring: **blocked**.
- Full LUMIERE four-analysis run: **not recommended**.
- Duke: **conditional future path; access reportedly granted, exact accepted terms not recorded, local archive incomplete, metadata inventory prepared and pending**.
- Fresh no-request cohort: **none verified in the bounded screen**.
- Two-history proposal: **r3 draft, not adopted; selector, shared-volume helper, and typed failures are synthetic-only; no patient runner is wired**.
- No project code initiated a download or external contact. The user reports a local download in progress. No archive member was opened, no patient voxel was read, and no scoring occurred. This plan does not authorize MRI/mask reads or patient scoring.

## Sources checked

- [Duke official record](https://zenodo.org/records/20548528)
- [RHUH-GBM official record](https://www.cancerimagingarchive.net/collection/rhuh-gbm/)
- [LUMIERE dataset record](https://springernature.figshare.com/articles/dataset/LUMIERE_dataset_-_MRI_data_and_automated_segmentations/21249516)
- [PREDICT-GBM paper](https://arxiv.org/html/2509.13360)
- [SoC-DT paper](https://arxiv.org/html/2510.03287)
- [UCSF-ALPTDG paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11294954/)
