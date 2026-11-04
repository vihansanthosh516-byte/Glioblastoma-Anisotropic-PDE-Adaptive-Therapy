# Track A and Track B full research plan

**Updated:** 2026-11-03
**Scope:** Tracks A and B only. Track C is outside this plan.
**Status:** Track A is closed for its current scope. Track B has a conditional Duke path. Metadata inventory is the only authorized data step. This document is a research plan, not a protocol amendment or authorization for MRI/mask reads or patient scoring.

## Executive decision

- **Track A:** close current work with A5 as an exploratory within-cohort result. It failed its frozen success rule. No screened outside cohort passed all requirements for independent confirmation. Do not keep scouting or run a new A analysis without a new question and user-approved scope.
- **Track B:** continue only through sequential Duke gates. Finish the archive download, verify integrity, and inventory approved metadata. Do not read MRI, mask, or landmark-coordinate contents or score patients at this stage.
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

### Proposed narrow question

For a patient with at least two prior scans, does PDE-based spatial ranking improve prediction of the next observed automated tumor-core mask over a two-history non-PDE shape forecast, when both methods receive the same history-only volume forecast, interval, and support?

This is a predictive comparison. It cannot isolate a causal effect of PDE biology. The target is agreement with the dataset's automated labels, not verified tumor truth. The horizon is the next observed visit conditional on its recorded interval, not a prospective fixed-day forecast.

This question and its comparator family remain **proposed**. No Duke scoring protocol is adopted. Do not copy LUMIERE's folds, timing, labels, cutoffs, or thresholds into Duke without a separate amendment.

### Minimum fairness requirements for a future Duke protocol

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

**Access and data status (2026-11-03):** The user reports download access and started downloading one `DATASET_FINAL.zip`. The local file was 76,144,640 bytes at the latest check, below the expected 17,043,569,844 bytes. The local Zenodo API snapshot gives MD5 `99580fb8ccb3daed186dc6afe5b5bf2f`. Download is incomplete. No ZIP member has been opened. No MRI, mask, landmark-coordinate, or patient-label voxel data has been read. No patient scoring has occurred.

The user authorized a metadata-only inventory covering filenames, anonymized patient/session identifiers, scan counts and dates, treatment/IDH metadata. This does not include MRI, mask, or landmark-coordinate contents. Exact accepted Duke license terms and written school approval are not saved. Zenodo lists CC BY-NC-ND 4.0 or a custom Duke agreement; the allowed use and sharing of derived masks, predictions, figures, and weights must be established from the accepted terms.

The prepared inventory is `src/fk_duke_archive_inventory.py`. It checks expected archive size and MD5, records SHA-256, lists central-directory metadata without opening members, and can read one explicit non-landmark CSV header. It writes locally under `output/track_b_duke_inventory_20261103/`; do not commit patient-level inventory outputs.

## Ordered work plan and gates

### Stage 1 — Complete and verify archive

1. Let the current download finish. Do not start a second copy or alter the partial file.
2. Check exact size and Zenodo MD5. Stop if either differs. Record local SHA-256 only after the expected checks pass.
3. List ZIP central-directory paths, sizes, and file types. Do not extract or open MRI, mask, or landmark members.

**Exit condition:** complete archive passes both expected size and MD5. A failed check stops the run.

### Stage 2 — Metadata-only cohort audit

1. Review the central-directory inventory and identify the dataset's metadata tables.
2. Inspect metadata headers first. Then read only approved patient/session IDs, relative scan dates, treatment and surgery dates, and IDH fields needed for eligibility and interval rules.
3. Do not read landmark coordinate tables, NIfTI files, brain masks, or segmentation masks in this stage.
4. Build a local, access-controlled patient/session map. Report aggregate counts and exclusion reasons; avoid publishing identifiers.
5. Determine the number of patients with usable three-visit sequences, valid intervals, eligible disease status, and no known project overlap. Do not call the cohort independent while overlap is unresolved.

**Exit condition:** the patient count, sample map, IDH status, treatment/surgery timing, interval policy, and overlap status are documented from allowed metadata. Set the minimum eligible cohort before any target-mask read. If it fails, stop the patient branch.

### Stage 3 — Rights and school approval

1. Save the written Duke access grant and exact accepted license/agreement.
2. Resolve whether analysis, derived masks, predictions, trained weights, figures, and summary results may be created or shared.
3. Save written school approval with date, approver, and conditions.

**Exit condition:** both written approvals and output-rights rules are recorded. Download access alone does not pass this gate.

### Stage 4 — Freeze Duke-specific protocol

Before any MRI or mask voxel read, write and adopt a dated Duke-specific amendment defining:

- target labels and endpoint, including treatment of cavity and empty masks;
- scan-pair/triplet selection and next-visit interval rule;
- origin-defined support and target-outside-support scoring;
- image registration method, QC thresholds, retry and exclusion rules;
- outer folds, inner selection, seeds, preprocessing, and patient-level estimand;
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

Track A's conclusion stays bounded: A5 completed within-cohort exploratory analysis and failed its frozen superiority rule. Track B, if completed, can support only a patient-level predictive comparison against source-provided automated masks in the eligible Duke cohort and observed intervals.

Neither track can claim a causal PDE effect, true tumor forecast accuracy, clinical benefit, or broad novelty from this design alone. Rigor, clarity, and reproducibility can improve through the gates above. Novelty and evidence depend on the results and the prior-art check. Do not assign a 7–8/10 score as a target or claim a positive result in advance.

## Key files

- Current Track A decision: `paper/trackA_novelty_decision_20261101.md`.
- Track B direction: `paper/track_b_codex_direction_and_execution_plan_2026-11-02.md`.
- Duke readiness: `paper/track_b_duke_readiness_2026-11-03.md`.
- Duke inventory code: `src/fk_duke_archive_inventory.py`.
- Historical LUMIERE go/no-go and control proposals remain historical or proposed; they do not authorize Duke scoring.
