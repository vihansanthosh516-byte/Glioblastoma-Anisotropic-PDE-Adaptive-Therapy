# Execution Plan for GBM Tracks A, B, and C

**Status:** Current project execution roadmap, recorded 2026-10-31. This file records the user's selected next work for all three tracks.

This roadmap adds execution order and stop rules. It does not replace a frozen protocol, change an endpoint, grant data access, or authorize patient-data scoring. Commit any required protocol amendment and code/input hashes before the next analysis. Keep raw data local and follow each source's terms.

## Shared rules

1. Keep Tracks A, B, and C separate. Their cohorts and endpoints differ. No patient-level link across tracks is established. Do not call the project a validated digital twin.
2. Label evidence as observed, exploratory, simulated, assumed, or post hoc. Keep negative results and state the weakest limitation beside each main result.
3. Treat novelty ratings as internal, subjective planning estimates. Current working estimates are A 3/10, B 3–4/10, and C 2–3/10.
4. Treat 6/10 as earned only when a focused prior-art review identifies a narrow gap and a frozen test supports it against strong controls. Treat 7/10 as a stretch requiring stronger independent confirmation or a clear new scientific insight. Do not promise either score or claim “first.”
5. A result that fails its frozen rule is still a completed result. Do not tune again after test outcomes are seen.

## Track A — donor-held-out molecular prediction and external signature test

**Current status.** The Yu 2020 donor-held-out run has no scores or `result.json`. The sealed A2 runner's timing pilot stopped before held-out prediction. Its conservative lower bound for all 14 planned outer folds is 10.89 hours, above the six-hour gate. A3 blocks scoring under the A2 grid. See `paper/trackA_yu2020/PROTOCOL.md` and `output/trackA_yu2020/pilot_timing.json`.

### A1. Complete the Yu 2020 analysis within the runtime gate

1. Prepare a dated protocol amendment and matching runner/checksum update before another expression-value read or any held-out scoring. First seek output-preserving runtime improvements. If these do not meet the gate, reduce the tuning grid or model scope in the amendment. Do not remove a model or change a setting silently.
2. Preserve donor-held-out splits, training-only filtering and scaling, the class-balanced AURC primary endpoint, the logistic-regression comparator, and the frozen descriptive success rule where feasible. The candidate qualifies as “better (descriptive)” only if mean AURC is lower and candidate AURC is lower for at least five of six primary donors.
3. Run a timing pilot for the exact amended scope. Include loading and overhead in the projection. Do not score held-out donors unless the projected full run fits the six-hour limit.
4. Target all eligible primary and sensitivity donors if the amended design passes the runtime gate. If it still fails, amend the scope before scoring: retain the six-donor primary cohort and label the eight-donor sensitivity cohort “not run.” Pilot the revised scope again.
5. Run once on the frozen design. Save `output/trackA_yu2020/result.json` only after all required folds finish. Read it back and verify hashes, donor list, cell flow, per-donor/per-cell-type metrics, runtime, and implementation choices. Report results unchanged, including a negative result.

### A2. Test a separate Core/Periphery signature for external transfer

This is a distinct study from the Yu 2020 classifier. Use the atlas only after confirming source labels, donor mapping, malignant-cell labels, assay compatibility, and reuse terms. Stop if any of these gates fail.

1. Finish a focused prior-art review. State a narrow question: does a donor-stable, malignant-cell Core/Periphery signature transfer to independent anatomical-region samples beyond normal-brain contamination and established cell-state signals?
2. Build the signature on atlas donors only. Use leave-one-donor-out checks; retain genes with the same direction in at least six of seven eligible donors, if seven pass source checks.
3. Before accessing Ivy GAP expression values, freeze the signature, patient-level endpoint, comparators, scoring code, and pass rule. Check data terms and obtain the user's approval before any download.
4. Score Ivy GAP once by patient. Compare with normal-brain contamination markers, Neftel cell-state scores, a prior Core/Periphery signature, and 1,000 size-matched random gene sets. Use patient-level bootstrap intervals. The candidate passes its stated rule only if it exceeds the contamination score with an interval above zero, exceeds the random-set 95th percentile, and remains associated after adjustment for contamination.
5. If the rule fails, report no external transfer. A pass could support a narrow novelty claim near 6/10; small patient count limits precision. A 7/10 claim needs independent replication or stronger evidence.

## Track B — patient-level PDE forecast with uncertainty and failure checks

**Current status.** No new patient-mask scoring is authorized. MU is previously examined and can support exploratory work only; source-backed case-level registration, label-zero meaning, and full-grid support evidence remain open. Duke is the strongest independent lead, but access and school-use approval remain pending. The preliminary B2.1 report and follow-up are committed as `44b25e5` and `0147183` in `paper/track_b_prior_art_2026-10-31.md` and `paper/track_b_prior_art_followup_2026-10-31.md`. SoC-DT already tests PDE plus neural correction on real UCSF patients, so that broad combination is not a novelty claim. The reports remain preliminary: direct PubMed/Scholar searches, citation sweeps, and several full-text/code checks remain open. The gap is unconfirmed.

### B1. Clear source and protocol gates before patient scoring

1. For Duke, wait for documented access approval and school/fair permission for the exact license and planned outputs. Do not download restricted files before both approvals.
2. For MU development, obtain source-backed case-level registration evidence, label-zero interpretation, and valid-support guidance. Keep any resulting patient scores exploratory.
3. Before opening target-mask voxels, create a metadata-only eligibility report. Check patient overlap and prior exposure, visit dates, mask provenance, registration quality, valid support, cohort flow, and achievable precision. Stop if a required source or quality gate fails.
4. Reconcile the endpoint before freezing the Track B protocol. The current research-plan appendix names future tumor-core geometry and patient-macro Dice as primary. The novelty refinement also proposes predicting the count and locations of newly annotated additions. Preserve future-core Dice as primary; use count/location scores as secondary diagnostics unless a dated amendment promotes them before outcomes are seen.

### B2. Test one specific candidate, not “PDE plus AI” generally

1. Close the remaining searches in `paper/track_b_prior_art_followup_2026-10-31.md`: direct PubMed/Scholar queries for glioma count models; citation sweeps for LUMIERE, MU-Glioma-Post, UCSF-ALPTDG, and PREDICT-GBM; full texts for the named PDE/residual methods and TumorFlow; remaining uncertainty and code/license checks; and any reported persistence Dice by horizon. Record what was checked and what remains unresolved. Until then, call the gap preliminary, not novel.
2. Freeze one primary target and horizon rule before scoring. Keep future tumor-core geometry and patient-macro Dice as primary under the current research plan. Treat new-addition count and location as secondary diagnostics unless a dated amendment promotes them before outcomes are seen. Report plain Dice at each predeclared forecast horizon. Never give any model the true future tumor size or use oracle volume matching for the primary score.
3. Keep the candidate low-dimensional: a tissue-aware PDE spatial prior, a training-only estimate of expected additions from prior scans, and a small residual correction shrunk when mask, registration, or tissue-map perturbations make the forecast unstable. Treat the stability score as a sensitivity measure, not calibrated biological uncertainty. Use the existing synthetic count/location gate only as feasibility evidence; verify the method at realistic region sizes before patient testing. Do not add an untested spatial-interaction model.
4. **Fixed-growth spec resolved 2026-10-31** (`paper/track_b_fixed_growth_baseline_spec_2026-10-31.md`; code `src/fk_fixed_growth.py`; 12 synthetic tests passed). One 76–104-day window; radius 0–20 mm in 2 mm steps, chosen on outer-training patients only, no inner folds; wired into the MU development runner; U-Net and count-model recipes frozen in spec §5 (not coded); distance-shell rule **dropped** (it was optional and never frozen, so no required baseline was removed); appendix §3.3/§3.5 amended. Valid support and registration remain source gates. Original requirement: require these mask baselines on the same held-out patients and horizon bins: (a) copying the last mask, (b) PDE-only, (c) a U-Net trained within development folds, and (d) fixed growth by dilation of the last mask. Tune dilation radius or growth size using inner training folds only, separately for each predeclared horizon bin. Construct each forecast without the target mask or future tumor size. The training-only count model remains in the comparison, but score it on count metrics only; a count alone has no location. The training-only distance-shell spatial rule was optional and is now dropped (2026-10-31). Also report the current plan's volume-trend, geometric-rate, and component-ablation comparisons. Do not use PREDICT-GBM's overlapping cohorts as fresh confirmation. Its fixed-volume coverage result is a warning, not a numeric threshold for plain Dice.
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

- **Track A:** draft and seal the runtime amendment; run the amended timing pilot before any held-out scoring.
- **Track B:** fixed growth is specified, coded, wired into the MU runner and synthetically tested; U-Net and count recipes are frozen but not coded. Next: code the U-Net and count model; source-backed valid support and MRI-only registration QC; cavity meaning from source label documents; finish Scholar and dataset citation sweeps. No patient-mask scoring until the source/access gates pass.
- **Track C:** finish the current negative/limited-results report; keep new model testing stopped until measured kinetics and independent review exist.

Track these actions separately. A pass in one track does not clear another track's data or protocol gates.
