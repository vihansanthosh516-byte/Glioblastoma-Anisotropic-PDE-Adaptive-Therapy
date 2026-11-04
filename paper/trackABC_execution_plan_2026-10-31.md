# Execution Plan for GBM Tracks A, B, and C

**Status:** Current project execution roadmap, updated 2026-11-04. Track B's current direction supersedes the earlier PDE-first and LUMIERE/MU plans below; those details are retained as historical record.

This roadmap adds execution order and stop rules. It does not replace a frozen protocol, change an endpoint, grant data access, or authorize patient-data scoring. Commit any required protocol amendment and code/input hashes before the next analysis. Keep raw data local and follow each source's terms.

## Shared rules

1. Keep Tracks A, B, and C separate. Their cohorts and endpoints differ. No patient-level link across tracks is established. Do not call the project a validated digital twin.
2. Label evidence as observed, exploratory, simulated, assumed, or post hoc. Keep negative results and state the weakest limitation beside each main result.
3. Treat novelty ratings as internal, subjective planning estimates. Current working estimates are A 3/10, B 3–4/10, and C 2–3/10.
4. Treat 6/10 as earned only when a focused prior-art review identifies a narrow gap and a frozen test supports it against strong controls. Treat 7/10 as a stretch requiring stronger independent confirmation or a clear new scientific insight. Do not promise either score or claim “first.”
5. A result that fails its frozen rule is still a completed result. Do not tune again after test outcomes are seen.

## Track A — donor-held-out molecular prediction; external confirmation unavailable

**Current status (2026-11-01).** The frozen A5 Yu 2020 analysis completed on its six primary donors: GS1, GS2, GS3, GS5, GS6, and GS13. The candidate's mean class-balanced AURC was lower than logistic regression (0.1479 vs. 0.1764), but it won in only four of six donors; the frozen success rule fails. Cell-type-only performance was within the frozen 0.05 margin of the best mapped-cell expression model. The eight-donor sensitivity analysis was not run. This is within-cohort donor-held-out validation, not independent external confirmation or a malignant-cell-only test. Results and limits: `output/trackA_yu2020/result.json` and `paper/trackA_execution_audit_20261031.md`.

The metadata-only screen found no new public, no-request cohort that clears all gates. SNUH has six pairs but prior project exposure; GSE286419 retains at most four source-validated regional pairs; GSE318564 has three explicit IDH-wildtype pairs and an unresolved annotation join; HRA009021 lacks a public sample-to-region crosswalk and verified reuse terms. Keep all outside data out of the current figure. Pause broad scouting until a specific new accession appears. Full comparison, novelty limits, controls, and stop rules: `paper/trackA_novelty_decision_20261101.md`.

### A1. Report the frozen result without retuning

1. Preserve the completed A5 result and protocol. Report the candidate's failed donor-level success rule, all required baselines and ablations, cell-type composition diagnostic, point-size diagnostic, and within-cell-type limits as recorded. Make no additional post-score tests, tuning, or biological claims.
2. Distinguish the six held-out donors inside the Yu cohort from an independent cohort. Do not describe A5 as malignant-cell-only validation or as external transfer.
3. Keep the internal novelty estimate unchanged. Current evidence does not support a higher claim. See `paper/trackA_novelty_decision_20261101.md`.

### A2. Reopen external validation only for a named lead

Require, before any download, scoring, or figure change: at least six IDH-wildtype patients with public metadata mapping both regions to the same donor; defensible cell-level malignant labels; clear no-request access and reuse terms; exact sample-ID checks against training and prior-use records; and verified independence from all previously used cohorts. Use donor as the unit, keep all samples/cells from each donor together, and freeze paired donor-level endpoints and controls before any expression read. Stop on any missing map, ambiguous label, unclear term, overlap, or donor count below six. Spatial spots alone do not satisfy the malignant-cell-label gate.

Do not use Ivy GAP, SNUH, or another exposed cohort as untouched confirmation. Do not resume broad scouting. A specific public accession may be screened metadata-only against every gate; if none appears, report external confirmation as unavailable.

## Track B — Duke longitudinal MRI biology and modeling (current; 2026-11-04)

**Current status.** The user reports Duke download access. The local archive passed the expected byte-size check and matched the MD5 supplied with the authorized download; a local SHA-256 was recorded. Metadata-only inventory is complete. No NIfTI member, MRI or segmentation-mask voxels, or landmark-coordinate table has been read, and no patient scoring has occurred. Exact accepted use terms, written school/competition review, prior-cohort overlap, image QC, and a frozen protocol remain open. See `paper/track_b_duke_readiness_2026-11-03.md`.

**Current proposed question.** In three consecutive visits, are middle-visit newly appearing nonenhancing FLAIR-labeled regions more likely than persistent regions at the same distance from the current core to receive a necrotic/enhancing-core label at the next visit? ADC prediction is secondary; reaction-diffusion forecasting is a separate secondary benchmark. This tests temporal transitions among automated MRI labels under observed care, not histologic infiltration, causation, clinical benefit, or a verified novelty claim. Full proposal: `paper/trackB_research_plan_draft_2026-11-03.md`.

### Current ordered work and gates

1. Save and review the written Duke grant and exact accepted terms for analysis, retention, and output sharing. Obtain written school/competition review.
2. Freeze the source-coded diagnosis/IDH eligibility rule, date-sorted window rule, minimum cohort and precision requirement; reconcile patient overlap with prior project cohorts where permitted.
3. Review and adopt a Duke-specific protocol. Keep the primary transition estimand separate from the secondary ADC and forecast comparisons. Freeze labels, support, registration, exclusions, folds, model settings, uncertainty, failure rules, and the runtime ceiling.
4. Complete synthetic integration and full-pipeline timing. If the measured projection exceeds the adopted ceiling or the cohort cannot support the frozen precision target, amend the scope before target-label scoring or stop.
5. After rights, school review, and a frozen protocol, perform the approved MRI-only registration/support QC with target segmentation labels closed. Freeze and hash the eligible cohort and QC list.
6. Seal code, settings, cohort, folds, and allowed input fingerprints before a single locked evaluation. Report all outcomes, missing arms, exclusions, and failures unchanged.

**Historical plans below.** The former PDE-first Duke proposal and earlier LUMIERE/MU work are not the current primary question and do not authorize any patient read or score. Their details remain below only as project history.

## Historical Track B — PDE forecast and LUMIERE/MU plans (superseded 2026-11-03)

**Historical status (2026-10-31 to 2026-11-02).** This paragraph records the former PDE-first question. Duke access is now reported by the user and the archive metadata inventory is complete; current gates and question are stated above. The older PDE prior-art screen remains preliminary and does not establish novelty for the current FLAIR-transition question.

### B1. Clear source and protocol gates before patient scoring

1. For Duke, wait for documented access approval and school/fair permission for the exact license and planned outputs. Do not download restricted files before both approvals.
2. For MU development, obtain source-backed case-level registration evidence, label-zero interpretation, and valid-support guidance. Keep any resulting patient scores exploratory.
3. Before opening target-mask voxels, create a metadata-only eligibility report. Check patient overlap and prior exposure, visit dates, mask provenance, registration quality, valid support, cohort flow, and achievable precision. Stop if a required source or quality gate fails.
4. Reconcile the endpoint before freezing the Track B protocol. The current research-plan appendix names future tumor-core geometry and patient-macro Dice as primary. The novelty refinement also proposes predicting the count and locations of newly annotated additions. Preserve future-core Dice as primary; use count/location scores as secondary diagnostics unless a dated amendment promotes them before outcomes are seen.

### B2. Test one specific candidate, not “PDE plus AI” generally

1. Close the remaining searches in `paper/track_b_prior_art_followup_2026-10-31.md`: direct PubMed/Scholar queries for glioma count models; citation sweeps for LUMIERE, MU-Glioma-Post, UCSF-ALPTDG, and PREDICT-GBM; full texts for the named PDE/residual methods and TumorFlow; remaining uncertainty and code/license checks; and any reported persistence Dice by horizon. Record what was checked and what remains unresolved. Until then, call the gap preliminary, not novel.
2. Freeze one primary target and horizon rule before scoring. Keep future tumor-core geometry and patient-macro Dice as primary under the current research plan. Treat new-addition count and location as secondary diagnostics unless a dated amendment promotes them before outcomes are seen. Report plain Dice at each predeclared forecast horizon. Never give any model the true future tumor size or use oracle volume matching for the primary score.
3. Keep the candidate low-dimensional: a tissue-aware PDE spatial prior, a training-only estimate of expected additions from prior scans, and a small residual correction shrunk when mask, registration, or tissue-map perturbations make the forecast unstable. Treat the stability score as a sensitivity measure, not calibrated biological uncertainty. Use the existing synthetic count/location gate only as feasibility evidence; verify the method at realistic region sizes before patient testing. Do not add an untested spatial-interaction model.
4. **Fixed-growth spec resolved 2026-10-31** (`paper/track_b_fixed_growth_baseline_spec_2026-10-31.md`; code `src/fk_fixed_growth.py`; 12 synthetic tests passed). One 76–104-day window; radius 0–20 mm in 2 mm steps, chosen on outer-training patients only, no inner folds; wired into the MU development runner; U-Net and count model coded per spec §5 (`src/fk_learned_baselines.py`) and wired into the same runner, 78 related synthetic tests passed; distance-shell rule **dropped** (it was optional and never frozen, so no required baseline was removed); appendix §3.3/§3.5 amended. Valid support and registration remain source gates. Original requirement: require these mask baselines on the same held-out patients and horizon bins: (a) copying the last mask, (b) PDE-only, (c) a U-Net trained within development folds, and (d) fixed growth by dilation of the last mask. Tune dilation radius or growth size using inner training folds only, separately for each predeclared horizon bin. Construct each forecast without the target mask or future tumor size. The training-only count model remains in the comparison, but score it on count metrics only; a count alone has no location. The training-only distance-shell spatial rule was optional and is now dropped (2026-10-31). Also report the current plan's volume-trend, geometric-rate, and component-ablation comparisons. Do not use PREDICT-GBM's overlapping cohorts as fresh confirmation. Its fixed-volume coverage result is a warning, not a numeric threshold for plain Dice.
5. Fit all tuning and calibration on training patients only. Report plain patient-macro Dice for each horizon, with paired patient intervals; report count calibration, location scores, surface distance, volume error, failures, and parameter identifiability as secondary outcomes. If a fair U-Net or fixed-growth baseline cannot be built, mark the novelty test incomplete; do not claim the candidate passed.
6. Freeze a meaningful pass margin and precision rule before external outcomes are opened. Test once on Duke only if every access and quality gate passes. If Duke fails or remains unavailable, keep MU exploratory and make no independent-confirmation or 6–7 novelty claim.
7. The candidate must beat every required mask baseline under the frozen rule. A tie with U-Net means the PDE prior adds interpretability, not demonstrated accuracy. Do not claim PDE plus learned correction as new; SoC-DT already reports that combination on real patients. Limit any remaining novelty claim to the training-only additions, stability-shrunk residual, and fair horizon-specific tests only if the remaining literature searches leave that gap open. Do not substitute CFB or UCSD annotation-change work for the PDE aim.

## Track C — close current simulation; reopen only with measured kinetics

**Current status.** Track C remains simulation-only, with assumed drug-response and recovery rates. Sealed H1 failed in tests 171 and 173. The planned 600-scenario primary outcome was not run. The 31 October decision stops further model testing. See `paper/trackc_final_status_2026-10-31.md`.

### C1. Complete the current study without changing its result

1. Report the current work as an exploratory mean-versus-worst-decile trade-off inside one assumed model family. Include both sealed H1 failures, the unrun primary outcome, the assumed-rate sensitivity results, and the limits of the simulator.
2. Do not change margins, endpoints, suites, or rate ranges. Do not run a new sealed test or revive the older headroom/probing simulation plan. The later stop decision supersedes that older proposal.
3. Do not claim measured drug response, patient benefit, or clinical treatment guidance.

### C2. Conditional future novelty path

1. Reopen model testing only if a source provides same-model off-drug recovery measurements at several prespecified washout durations, preferably in patient-derived GBM cells, with replicate-level post-washout TMZ dose-response.
2. Ask an independent reviewer to confirm that these measurements bound the switching and reversion rates before policy outcomes are examined. If no suitable evidence exists, keep Track C closed as a simulation study.
3. If the evidence gate passes, write and freeze a separate protocol. Use fresh scenario seeds, a power/precision justification for the tail outcome, and predeclared comparisons against the best fixed and simple adaptive schedules. Preserve all earlier failures; new evidence does not erase them.
4. Treat this as a calibrated simulation result, not patient benefit. A direct patient-benefit study needs an unexposed patient-level dataset linking actual TMZ dose and timing, decision-time tumor state, toxicity, and later outcomes. The current search found no qualifying open dataset. Do not claim a 6–7 novelty rating until an eligible source and independent test exist.

## Next actions

- **Track A (updated 2026-11-01 audit):** A5 is complete and failed its frozen 5/6 donor rule (`output/trackA_yu2020/result.json`). Report external confirmation as unavailable. Reopen only for a named accession that passes every gate (`paper/trackA_novelty_decision_20261101.md`). The earlier "seal the runtime amendment and run the timing pilot" action is done.
- **Track B (current Duke path):** save and verify the accepted data-use agreement and written school/competition review; resolve permitted patient overlap and metadata eligibility; review and freeze the biology-first protocol; complete synthetic integration and runtime timing; then, only after rights and protocol gates, conduct approved MRI-only QC and freeze the cohort before target-label scoring. Current status: `paper/track_b_duke_readiness_2026-11-03.md`.
- **Track B historical LUMIERE path (2026-11-01; superseded for the current question):** the exploratory plan and Amendment R are documented in `paper/track_b_no_contact_dataset_plan_2026-10-31.md`; no LUMIERE scoring occurred. Its 677.6 wall-clock-hour projection and open gates do not describe the current Duke path.
- **Track B historical MU path (2026-10-31):** fixed-growth, U-Net, and count code and the historical runtime/gate record remain in `paper/track_b_runtime_pilot_2026-10-31.md` and `paper/track_b_amendment_decision_2026-10-31.md`. These MU settings and projections are not the current Duke protocol.
- **Track C:** finish the current negative/limited-results report; keep new model testing stopped until measured kinetics and independent review exist.

Track these actions separately. A pass in one track does not clear another track's data or protocol gates.
