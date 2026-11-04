# Track A and Track B full research plan

**Updated:** 2026-11-04
**Scope:** Tracks A and B only. Track C is outside this plan.
**Status:** Track A is closed for its current scope. Track B is being reframed around a biology-first Duke study. The downloaded archive passed its expected size and MD5 checks and a metadata-only inventory was completed. No NIfTI member or landmark-coordinate table was opened; no MRI or mask voxel array was read; no patient scoring occurred. This document is a research plan, not an adopted protocol or scoring authorization.

## Executive decision

- **Track A:** close current work with A5 as an exploratory within-cohort result. It failed its frozen success rule. No screened outside cohort passed all requirements for independent confirmation. Do not keep scouting or run a new A analysis without a new question and user-approved scope.
- **Track B:** make the three-visit Duke tissue-state transition question primary: compare next-visit core-label conversion in newly appearing versus persistent FLAIR regions at matched distance. Treat ADC as secondary and PDE forecasting as a secondary benchmark. Finish the research plan and protocol freeze before image/mask scoring; exact accepted data-use terms and school review remain open.
- **No positivity target:** A null or negative result can be a valid study. Do not tune or choose research questions to obtain a positive score. Novelty and evidence cannot be promised as 7–8/10 before the required data and independent evaluation exist.

## Track A: closeout

### Question and result

A5 tests whether the candidate ranks confidently classified cells better than its comparator under donor-held-out validation. Its target is the project-defined MRI-marked region label in mixed-cell expression data. The six donors, not the thousands of cells, are the independent units.

A5's candidate had lower mean AURC than logistic regression but won in only four of six donors. The frozen rule required at least five donor wins. Therefore A5 did not demonstrate superiority. Its result remains exploratory and within-cohort. It does not establish malignant-cell-specific biology, external validity, mechanism, or clinical value.

### Track A decision

Track A is **done for the current scope**. Preserve the frozen A5 analysis and report. Do not run the deferred calibration, within-library, scANVI, or attention analyses unless a dated scope amendment authorizes them. Do not resume broad dataset scouting.

Reopen Track A only if the user chooses a new question or a candidate cohort clears its stated gates. Label any smaller or reused cohort exploratory. Call a cohort independent only after patient-level overlap is resolved. The previous six-donor threshold applies to the project's external-confirmation claim; a smaller cohort would need a different, explicitly exploratory question.

### Track A records

- Frozen protocol and result: `paper/trackA_yu2020/PROTOCOL.md` and `output/trackA_yu2020/result.json`.
- Dataset screen and stopped leads: `paper/trackA_alternate_dataset_plan_20261101.md`.
- Novelty and closure decision: `paper/trackA_novelty_decision_20261101.md`.
- Do not change A5 results or figures to improve the narrative.

## Track B: research question

### Primary biology question

Among patients whose released records identify glioblastoma and code IDH mutation status as Negative, are newly appearing nonenhancing FLAIR-labeled regions more likely than persistent regions at the same physical distance from the current core to receive a necrotic- or enhancing-core label at the next visit? Does local ADC trajectory add predictive information?

The primary outcome is a distance-standardized, patient-level difference in automated-label conversion. The ADC analysis is secondary. A patient-specific reaction-diffusion forecast and its controls are a separate secondary benchmark; forecast performance will not decide whether the primary transition question is supported. This is a radiographic tissue-state study, not a histology test or a claim that FLAIR marks malignant cells.

Prior work already maps future recurrence-prone versus stable peri-enhancing tissue and reports ADC differences, and serial FLAIR and ADC progression analyses exist. The candidate distinction is the three-visit classification of middle-visit FLAIR as newly appearing or persistent, compared at matched distance with patient-level next-visit core-label conversion. This remains a candidate extension, not a verified novelty claim. See the competition-form draft and its references for the complete rationale.

This biology-first question and the detailed methods remain **proposed**. No Duke scoring protocol is adopted.

### Fairness requirements for the secondary forecast benchmark

- Use patient-grouped outer evaluation. Keep every patient and all visits in one fold.
- Fit volume forecasts, preprocessing, and comparator settings on outer-training patients only. Nested selections stay inside each outer-training set.
- Give every arm the same eligible patients, forecast interval, origin support, and history-only volume draws where those quantities apply.
- Compare the PDE candidate with frozen simple baselines and the proposed two-history control. Include PDE-only and any learned baseline only if the Duke-specific protocol defines them, their recipes are reproducible, and runtime is feasible.
- Define exclusions and recoverable failures before scoring. Report all fallbacks and keep every eligible patient in the denominator unless a frozen rule says otherwise.
- Freeze one primary metric, patient-level aggregation, uncertainty method, multiplicity family, pass rule, precision requirement, and compute ceiling before target labels are opened.
- Treat target masks outside origin-defined support as misses. Never derive support or tune registration from target labels.

The earlier synthetic two-history screen is a diagnostic only. It did not establish a strong patient baseline and cannot choose the winning settings. The revised selector must remain fold-local; no stability rule may use outcomes across outer folds.

## Track B: Duke evidence and current state

The [Zenodo record](https://zenodo.org/records/20548528) describes 61 patients and 503 sessions, with 3–18 sessions per patient. It lists atlas-standardized 1 mm MRI and automated segmentations. Its labels are 1 necrotic tumor core, 2 nonenhancing T2-FLAIR abnormality, 3 enhancing tumor, and 4 resection cavity. These cohort totals do not establish eligible independent triplets.

The [author usage guide](https://github.com/deshanyang/GBM-DIR-QA) says scans are atlas-aligned but separate sessions are not registered to each other. It also says dates were automatically extracted and may be offset by days. A longitudinal study therefore needs a frozen per-patient registration and interval policy. Automated masks remain source labels, not ground truth.

**Access and data status (metadata audit completed 2026-11-03; document updated 2026-11-04):** The user reports Duke download access. The local `DATASET_FINAL.zip` reached 17,043,569,844 bytes and matched the MD5 supplied with the authorized download (`99580fb8ccb3daed186dc6afe5b5bf2f`); a local SHA-256 was recorded. The restricted Zenodo API did not expose a file checksum for independent comparison. The ZIP central directory and permitted patient/session, relative-date, diagnosis/IDH, and treatment metadata were inspected. No NIfTI member, MRI or mask voxel array, or landmark-coordinate table was opened. The inventory remains outside the repository. No patient scoring occurred.

The metadata-only inventory stage is complete. It covered filenames, anonymized patient/session joins, relative scan dates, diagnosis/IDH metadata, and treatment-event metadata; it did not include MRI, mask, or landmark-coordinate contents. Exact accepted Duke license terms and written school approval are not saved. The draft proposes to include only records coded as glioblastoma with IDH `Negative`, while excluding missing, indeterminate, or mixed diagnoses; this is not comprehensive IDH1/2 confirmation. Prior-use overlap is not reconciled, so independence is unestablished. Zenodo lists CC BY-NC-ND 4.0 or a custom Duke agreement; allowed analysis and sharing of derived masks, predictions, figures, and weights must be confirmed from the accepted terms.

The metadata-only inventory tool is `src/fk_duke_archive_inventory.py`. It checks expected archive size and the MD5 supplied with the authorized download, records SHA-256, lists central-directory metadata without opening members, and can read one explicitly named non-landmark CSV header. The inventory record is saved outside the repository; do not commit patient-level inventory contents.

## Ordered work plan and gates

### Stage 1 — Archive verification (complete)

The completed ZIP matched its expected byte size and the supplied MD5. The inventory tool also recorded a local SHA-256 and listed only ZIP central-directory metadata. The archive passed this gate. No file member was extracted or decompressed.

### Stage 2 — Metadata-only cohort audit

The authorized metadata tables were read and session keys were matched to archive paths. The draft now states how source-coded IDH and mixed diagnoses will be handled; those rules and the resulting eligible patient list have not been frozen. Prior-use overlap is not reconciled, and image QC has not been run. Keep patient-level metadata local and do not publish it.

**Exit condition:** before image/mask reads, confirm and freeze the proposed IDH/diagnosis interpretation, date-sorted window rule, prior-use overlap status, and minimum-cohort/precision rule. If eligibility, independence, or rights fail, stop the patient branch.

### Stage 3 — Rights and school approval

1. Save the written Duke access grant and exact accepted license/agreement.
2. Resolve whether analysis, derived masks, predictions, trained weights, figures, and summary results may be created or shared.
3. Save written school approval with date, approver, and conditions.

**Exit condition:** both written approvals and output-rights rules are recorded. Download access alone does not pass this gate.

### Stage 4 — Freeze Duke-specific protocol

Before any MRI or mask voxel read, write and adopt a dated Duke-specific amendment defining:

- primary endpoint: next-visit conversion of newly appearing versus persistent nonenhancing FLAIR at matched distance; operational labels, cavity, and empty-mask rules;
- scan-pair/triplet selection and next-visit interval rule;
- origin-defined support and target-outside-support scoring;
- image registration method, QC thresholds, retry and exclusion rules;
- the primary patient-level transition estimand and bootstrap; separate outer folds, inner selection, seeds, and preprocessing for the secondary forecast and ADC prediction analyses;
- candidate and every required baseline, including the two-history control if retained;
- expected failures, common fallback, denominator, missing-arm rules, and incomplete outcomes;
- primary metric, secondary metrics, bootstrap/interval method, multiplicity, pass rule, and precision requirement;
- runtime ceiling, checkpointing/resume behavior, and plan if the estimate exceeds the ceiling.

Do not use the prior LUMIERE amendment as the Duke protocol. Freeze code, settings, cohort rules, and the primary contrast before outcomes.

**Exit condition:** user adopts the amendment. Any unresolved item keeps patient scoring blocked.

### Stage 5 — Image-only preflight and support validation

Only after a separate explicit read scope:

1. Load only the MRI inputs required for the frozen registration/support QC.
2. Register each longitudinal visit with the frozen MRI-only method. Keep target segmentation masks closed.
3. Validate geometry, transforms, support, finite values, and predeclared QC criteria. Use landmark pairs only if their use is separately permitted and frozen in the protocol.
4. Save QC metrics, file hashes, pass/fail reasons, and exclusions. Save no voxel arrays or images unless separately allowed.
5. Freeze and hash the eligible patient list before any target-mask access.

**Exit condition:** registration/support checks pass for enough patients to meet the frozen minimum. If not, stop; do not relax thresholds after seeing failures.

### Stage 6 — Synthetic integration and runtime gate

Build a Duke-specific loader and runner only after Stages 2–4 fix its schema and protocol. Keep the existing LUMIERE runner separate. Test synthetic data for fold isolation, target-read guards, common patients/support, history-only fitting, fallback accounting, and result hashing.

Then measure the complete Duke pipeline on synthetic data using the actual grid, interval distribution, folds, and all required arms. Count calibration, nested selection, U-Net training, refits, prior analyses, and reruns. Use the measured estimate and uncertainty to compare with the compute ceiling. LUMIERE wall-clock projections are not Duke estimates.

**Exit condition:** all required arms run end-to-end on synthetic data, and measured runtime fits the adopted budget. Otherwise amend scope before any patient outcome is inspected or stop.

### Stage 7 — Seal and locked evaluation

After all prior gates pass, create a pre-run seal at a clean commit. Include code and protocol hashes, package/runtime versions, cohort hash, settings, folds, input hashes allowed by the read scope, and authorization record. Confirm the sealed commit matches the current commit.

Run the frozen evaluation once. Do not tune on test outcomes. Report all required arms, patient-level paired differences and intervals, exclusions, failures, fallbacks, runtime, and negative findings. If any required arm is missing, label the result incomplete. Do not turn an incomplete or negative result into a pass.

## Stop rules

Stop the Duke patient branch if the archive fails integrity checks; written access, output rights, or school approval fail; overlap remains unresolved; too few independent eligible patients remain; the frozen registration/support QC fails; the required arms cannot be run fairly; runtime exceeds the adopted ceiling; or code, settings, cohort, and protocol hashes do not match.

No fallback to LUMIERE as "independent confirmation." No post-outcome threshold changes. No positive-result promise.

## Claims and quality targets

Track A's conclusion stays bounded: A5 completed within-cohort exploratory analysis and failed its frozen superiority rule. Track B, if completed, can support a patient-level association between history-defined FLAIR states and next-visit automated core labels, a separate held-out ADC prediction result, and a distinct forecast-benchmark result in the eligible Duke cohort under observed intervals.

Neither track can claim a causal PDE effect, true tumor forecast accuracy, clinical benefit, or broad novelty from this design alone. Rigor, clarity, and reproducibility can improve through the gates above. Novelty and evidence depend on the results and the prior-art check. Do not assign a 7–8/10 score as a target or claim a positive result in advance.

## Key files

- Current Track A decision: `paper/trackA_novelty_decision_20261101.md`.
- Track B current competition-form science plan: `paper/trackB_research_plan_draft_2026-11-03.md`.
- Track B current readiness and gates: `paper/track_b_duke_readiness_2026-11-03.md`.
- Historical PDE-first direction memo: `paper/track_b_codex_direction_and_execution_plan_2026-11-02.md`.
- Duke inventory code: `src/fk_duke_archive_inventory.py`.
- Historical LUMIERE go/no-go and control proposals remain historical or proposed; they do not authorize Duke scoring.
