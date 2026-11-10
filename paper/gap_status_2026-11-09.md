# Project evidence table, action plan, and claim boundaries (2026-11-09, second pass)

This record replaces the first-pass gap list of the same date. It was built from four independent read-only reviews (Track A cohorts, Track B records and terms, Track C settings drift, integrated claims), each checked against the files, plus one predeclared synthetic screen and one user-authorized runtime timing. No frozen rule, endpoint, threshold, cohort, or saved result was changed. No Track C policy outcome was generated. No Duke endpoint score was computed.

Status words: **fixed** (closed in this pass), **recorded** (the files already state it correctly), **blocked** (needs a user document, a user decision, or work that must not start yet), **irreducible** (no code change removes it).

## 1. Evidence table

### Whole project and integration

| # | Gap | Status | Source | Root cause | Fixable? |
|---|---|---|---|---|---|
| W1 | No matched patient-level biology, MRI, and drug response | Irreducible | `tracks/integrated/README.md` "How the tracks connect"; `configs/trackABC_integrated_v2.json` `input_rules` | The three datasets are different patients; none measures drug response | No |
| W2 | Tracks A and B inform Track C only as motivation | Recorded | Master plan line 126 (hypothetical molecular-like and MRI-like observations); line 56 and config `input_rules` (observed A5 metrics are not inputs); line 127 (`a, b, d, q` are simulation assumptions) | Design choice forced by W1 | No. The only observed-data channel is an optional Duke history-only aggregate for the MRI scale; it is unapproved, so the standard-normal fallback applies |
| W3 | Active master plan is stale on recorded outcomes | Fixed (status note added; plan text untouched) | Master plan lines 19-21, 25, 88, 143 describe A5 and the 119-gene proxy as future; line 37 gives a directional MGMT hypothesis | Plan written in future tense for submission and never updated | Yes, by the note. Choosing the revised 2026-11-08 plan as the active record is a user decision (U5) |
| W4 | Revised 2026-11-08 plan omits the Student-t failure and the Track C stop | Fixed (status note added) | Revised plan lines 16, 265 keep the 14-day positive rule as if attainable | Written before the feasibility audit | Yes |
| W5 | Missing success/failure rules for analyses already run | Irreducible | Track B: master line 31 gives a direction only. Track C observational: line 37 direction only | No threshold was declared before the outcomes were seen | No. A rule cannot honestly be added after the result is known. These analyses are reported as descriptive |
| W6 | Title says "therapy optimization" | Recorded as a claim limit | Master plan line 5 | Title kept by the user | Wording choice for the user; every summary must say "in-silico" and "simulated" |

### Track A

| # | Gap | Status | Source | Root cause | Fixable? |
|---|---|---|---|---|---|
| A1 | A5 failed its frozen 5-of-6 donor rule | Irreducible; preserved | `output/trackA_yu2020/result.json`; `START_HERE.md` Yu 2020 bullet: candidate mean AURC 0.1479 vs logistic 0.1764; wins in 4 of 6 donors | The candidate did not win consistently across donors | No |
| A2 | No gain over the no-gate version; worse Brier | Irreducible | Same bullet: no-gate AURC 0.1473 vs 0.1479; Brier 0.2335 vs 0.1944 | Same | No |
| A3 | 119-of-121 gene score is a published-response proxy | Irreducible | `paper/trackA_tmz_response_119of121_amendment_2026-11-08.md`; paired differences split 3 and 3 across the six donors | In-vitro signature applied to mixed-cell region labels | No |
| A4 | Labels are sampled regions, not verified malignant cells | Irreducible | `tracks/track-a/README.md` interpretation boundary | Source data design | No |
| A5 | No independent cohort passes every gate | Blocked (no dataset exists) | Gates: `paper/trackA_alternate_dataset_plan_20261101.md` lines 5, 75-80; `paper/trackA_dataset_audit_20261031.md` lines 59, 62, 115. New search 2026-11-09, below | Public multi-region GBM single-cell cohorts are small | Not by this project without contacting data owners |
| A6 | Calibration, within-library, scANVI, attention analyses | Blocked; stay deferred | `START_HERE.md` Track A status bullet | Need a dated amendment before any run | Only by a user-approved amendment |

Track A cohort search, 2026-11-09 (metadata pages only; nothing downloaded). No candidate passes. Thirty-plus earlier screens are in the alternate-dataset plan. New candidates checked:

| Candidate | Verified facts | Gate failed |
|---|---|---|
| Wang et al. 2024, Sci Adv (CNGB CNP0003766 / NGDC PRJCA025943) | 5 patients, paired core and peritumoral brain, snRNA plus snATAC, malignant cells annotated ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11584018/)) | 5 donors, below 6; per-donor IDH status and access terms not verified |
| Moniz Garcia et al. 2026, Neuro-Oncology | 6 IDH-wildtype cases, MRI-guided core and margin samples ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13550674/)) | Bulk RNA, no malignant-cell labels; data on request |
| GSE253080 | 16 Visium samples from 3 patients ([GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253080)) | 3 donors |
| GSE242838 | 23 bulk samples by PET metabolism | Not a core/periphery map |
| GSE138794 | 28 gliomas; region inferred from Ivy GAP | Region label is expression-derived; raw data not available |

Nearest miss remains Pai et al. 2025 (GSE286419 / GSE286415): 4 source-validated Core/Margin pairs, one IDH1-mutant tumor. Not searched in this pass: Synapse, EGA public metadata, ArrayExpress. Absence of a candidate is therefore not proven.

### Track B

| # | Gap | Status | Source | Root cause | Fixable? |
|---|---|---|---|---|---|
| B1 | Duke cohort already scored | Recorded and verified | Private run records (outside Git): protocol freeze 2026-11-04 08:26 UTC; target-read gate 13:49; result 13:56 with status `EXPLORATORY_COMPLETED`, `target_voxels_read=true`, `confirmatory_threshold_met=false`; 33 patients, 101 windows | The primary endpoint was run once on the whole eligible cohort | No. The cohort is spent |
| B2 | Further exposure of the same cohort | Recorded | A second private analysis (label destinations) also read target masks; a third (ADC increment) has no result file and an incomplete preflight | Post hoc exploration | No. All of it is exploratory |
| B3 | Repo 2-mm code is not the scored analysis | Recorded and verified | Scored run: persistent = within 10 mm of previous FLAIR; excludes voxels within 10 mm of core; three bands, all beyond 10 mm. `src/track_b_transition_estimand.py`: 2-mm persistence (line 32), no 10-mm core exclusion (lines 28-30), four bands starting at 0 mm (line 37), different band-edge handling | Two different drafts | Not a defect to fix. The 2-mm code is a separate estimand. Running it on this cohort would be exploratory and is not planned |
| B4 | Outside-field voxels counted as non-converted | Irreducible for the old result | `src/fk_longitudinal_prep.py` line 100 (`apply_rigid`, constant fill 0); no coverage mask in the scored script | Resampling choice in the scored run | The old result keeps this caveat. The repo prototype has a coverage helper for any future, different cohort |
| B5 | Cohort contents | Recorded | Zenodo record (live, 2026-11-09): 61 patients, 503 sessions, 3 to 18 sessions each, six modalities, labels 1-4; archive listing: 3,521 NIfTI, 246 CSV, 1 workbook | - | Done |
| B6 | Written data-use terms | Blocked on the user (U1-U3) | No licence, readme, or terms file exists inside the archive (member-name scan). The live Zenodo page showed "Restricted" and a request form; the licence text was not visible to an unauthenticated read. The GitHub guide has no licence text | The accepted agreement was never saved | Only by the user saving the document |
| B7 | Registration, support, label, QC rules frozen before scoring | Blocked | No untouched Duke cohort exists for the primary endpoint | B1 | Only with a fresh cohort |
| B8 | Two-history control is not a strong comparator | Recorded | `output/track_b_two_history_screen_20261102/screen.json`: default control 0.8854 vs simpler control 0.8924 (growth), 0.7607 vs 0.7965 (shrinkage); synthetic, 27 cases per regime | The control does not beat last-mask extrapolation | No; it is labelled |
| B9 | Runtime for a Duke case | Fixed (timed) | One user-authorized case, four history files, no target file, no score: one rigid registration 13.0 s, label resample 0.24 s, state construction 0.20 s, total 13.6 s on a 182 x 218 x 182 grid (private record, outside Git) | - | Done. Two registrations per window gives about 27 s; 219 candidate windows would take about 1.7 hours. One case only; variance unknown |
| B10 | Fresh longitudinal cohort | Blocked (none found) | Search 2026-11-09, below | Public post-treatment cohorts lack three or more labelled visits, or were already used | Not at present |

Track B fresh-cohort search, 2026-11-09. Already used or scored by this project, so not fresh: Duke, MU-Glioma-Post, LUMIERE, RHUH-GBM, TUM-GBM, PREDICT-GBM, UCSD-PTGBM, CFB-GBM, BraTS material, UCSF. Candidates:

| Candidate | Verified facts | Why it does not qualify |
|---|---|---|
| Burdenko-GBM-Progression ([TCIA](https://www.cancerimagingarchive.net/collection/burdenko-gbm-progression/)) | 180 patients, 1 to 8 follow-ups, 131 GB | Only radiotherapy target volumes for a subset; no FLAIR or enhancing labels; controlled-access policy |
| BraTS 2024 post-treatment (Synapse syn59059776) | Catalogue listing only | Visits per patient unverified; registration required; BraTS material already examined |
| OpenNeuro ds007045 | 363 primary GBM patients, 17.8 GB, automatic masks | No longitudinal structure stated |
| Brain-Tumor-Progression (TCIA) | Catalogue listing: 20 patients | Too small; visits unverified |

Compute would be feasible for a cohort of Duke's size (B9). Data are the limit, not compute. Not checked: BraTS 2025 post-treatment, GLIS-RT, ReMIND, Yale collections.

### Track C

| # | Gap | Status | Source | Root cause | Fixable? |
|---|---|---|---|---|---|
| C1 | Assumed biology and signal quality | Irreducible | Master plan lines 121, 127 | Simulation | No |
| C2 | Progression and stop rules are simulation definitions | Irreducible | Config `stop_rule_status`: "synthetic proxy, not formal RANO" | Simulation | No |
| C3 | No validated interval method | Partly fixed | [GE-betting audit and screen result](trackc_ge_betting_coverage_audit_result_2026-11-09.md): code matches the published equations; GE-betting passed 30/30 synthetic cells (joint coverage 95.95% to 100%); Student-t failed 18/30 | Earlier methods were approximate | A valid candidate now exists. Still needed: independent statistical review and a materialized independent world schedule. Not adopted |
| C4 | 14-day precision target | Irreducible at any feasible world count; target kept | Same result: 0 of 30 cells reached 14 days for any method. GE-betting median half-width 18.3 days at 207 worlds per family with SD 30 days; 400.5 days at 7 worlds | Bounded-mean precision needs hundreds to thousands of worlds; the seed block allows 7 per family | Not without hundreds of runtime hours. Confirmatory evaluation stays stopped |
| C5 | "207 worlds" misread as sufficient | Fixed (documented) | Same result, section 3 | 207 only clears one worst-case lower bound | Done |
| C6 | Runtime | Recorded with stated uncertainty | One timing observation: 0.793 hours per world pair at 2 patients per cohort; 164 hours for 207 worlds per family | Single pilot | Error is unknown; cohort size multiplies it |
| C7 | Seed schedule, driver, sensitivity and bootstrap orchestration | Blocked by C4 | Config `not_yet_locked_before_evaluation`; `src/trackc_world_runner.py` is a development prototype with no entry point; only tests call it | Not built | Building them is pointless while C4 stops evaluation |
| C8 | Config is not loaded by any module | Partly fixed | Review: no source file reads config keys at run time; `efc4ea6` compares sampler axes, q levels, patient axes, sigma levels, and the stop code | Values are duplicated in code | New drift test `tests/test_trackc_config_parity_extended.py` now pins 9 further groups: clock, scan interval, particle-filter prior, resistance-score reference, controller settings, outcome settings, precision target and level, difference support, rejected-method status. A runtime loader still belongs in the driver |
| C9 | Controller model grid disagrees with config | Blocked on a user decision (U7) | `src/212_trackc_signal_fisher_grid.py` line 30: 3 q levels, 81 models; config `synthetic_q_signal_grid` has 5 levels, `candidate_model_set_maximum` 135 | Config grid was widened after the code | Flagged by a strict expected-failure test. Which grid is intended is a design choice, not made here |
| C10 | 112-day scan sensitivity cannot run | Blocked on a user decision (U7) | `PolicyBenchmarkSettings` rejects `scan_every=4` with 6 forecast cycles (`src/211_trackc_policy_benchmark.py` line 86) | Sensitivity declared without a forecast length | Flagged by a strict expected-failure test |
| C11 | Direct `RobustBeliefMPC` default stop rule is the legacy rule | Recorded | `src/208_trackc_robust_belief_mpc.py` line 117 default `"either"`; the benchmark overrides it | Historical default kept so old outputs reproduce | Flagged by a strict expected-failure test; the driver must pass the rule explicitly |
| C12 | Settings with no code counterpart | Blocked by C7 | Scan missingness 0.0 / 0.1, course-penalty sweep 0 / 10, bootstrap resamples, out-of-grid anchors of 0.75 | Orchestration not built | With the driver |
| C13 | Rejected max-T helper still callable | Recorded | `src/trackc_world_summary.py` line 171 | Kept for its audit script and tests | The parity test pins its "rejected" status in the config; it must not be used for analysis |

### Housekeeping

| # | Gap | Status | Source | Fixable? |
|---|---|---|---|---|
| H1 | Dirty worktree | Recorded | `git status`: 14 modified tracked files in `output/patient_split/` and one prior-art file, plus many untracked folders. None touched or staged | Preserved |
| H2 | Zenodo token pasted in chat | Blocked on the user (U4) | No token found in tracked or searched untracked text | Only the user can revoke it |

## 2. Ranked action plan

**Done in this pass (no rule changed):**

1. Audited GE-betting against the paper and ran a predeclared outcome-free coverage screen against three other methods.
2. Added a settings-drift test for values `efc4ea6` did not cover; three known drifts are pinned as strict expected failures so they cannot pass unnoticed.
3. Verified the Duke prior-use record, the 10-mm rule, and every difference from the 2-mm code.
4. Timed one Duke case under the user's authorization, without target files or scores.
5. Searched for fresh Track A and Track B cohorts.
6. Added status notes to both plans and updated the plan index.

**Needs the user (documents or decisions):**

- U1. Save the accepted Duke licence option and its written terms in `paper/`.
- U2. Save the school approval record.
- U3. Decide, from that text, whether Duke-derived numbers, figures, or model weights may be shared. Reported so far: the user states the file permits adapted versions and sharing; neither statement is verified from a document.
- U4. Revoke the pasted Zenodo token.
- U5. **Decided 2026-11-09:** the user made the revised 2026-11-08 plan the active plan. The 2026-11-07 master plan is superseded; the submission copy is form text only.
- U6. Track C: either keep confirmatory evaluation stopped (recommended) or authorize a fixed world count with a runtime cap, knowing the run would very probably end "inconclusive" under the unchanged 14-day rule.
- U7. Track C design choices that code cannot settle: 81 versus 135 controller models, and the forecast length for the 112-day scan sensitivity.

**Limits no code change removes:** W1, W5, A1-A4, B1-B4, C1, C2, C4.

## 3. Each track's question, data, procedure, endpoint, comparison, unit, rule, and claim boundary

| | Track A | Track B | Track C (simulation) | Integrated |
|---|---|---|---|---|
| Question | Does the A5 candidate improve donor-held-out Core-versus-Periphery prediction over L2 logistic regression? | Do newly appearing FLAIR voxels convert to core at a different rate than persistent FLAIR voxels at matched distance? | In synthetic worlds, does a dependence-aware minimax controller have higher mean utility than a matched `d = 0` controller? | Same as Track C, with molecular-like and MRI-like signals whose quality and dependence are assumed |
| Data | GSE117891, six strict donors, mixed-cell expression | Duke Zenodo 20548528, 33 patients, 101 windows (private) | Synthetic worlds only | Synthetic only; no patient values |
| Procedure | Frozen nested donor-held-out CV (`paper/trackA_yu2020/PROTOCOL.md`) | Private frozen protocol of 2026-11-04, 10-mm rule | Draft world benchmark; not frozen | Config v2; not frozen |
| Endpoint | Class-balanced AURC | Patient-averaged difference in conversion fraction, three bands | Utility: min of two progression times capped at 730 days, minus 5 days per course | Same |
| Comparison | L2 logistic regression (primary); kNN, SVM, GBT | New versus persistent FLAIR in the same patient and band | Matched `d = 0` controller; null-information controls | Same |
| Validation unit | Donor | Patient | Independent world | Independent world |
| Success / failure rule | Lower mean AURC and wins in at least 5 of 6 donors | None declared beyond the private protocol's confirmatory threshold | Both familywise 95% lower bounds above zero and both half-widths at most 14 days; otherwise feasibility or inconclusive | Same |
| Outcome | **Failed**: 4 of 6 | Exploratory, completed; private record states the confirmatory threshold was not met | **Not run; stopped** (precision infeasible) | **Not run** |
| Claim boundary | Region-label prediction in mixed cells within one cohort. No resistance, malignant-cell, or external claim | Exploratory description on one exposed cohort with a coverage caveat. No fresh confirmation | Simulated outcomes under assumed biology. No clinical or patient claim | A simulation design and a feasibility finding. Not a validated therapy framework |

## 4. Strongest defensible integrated question

"Across independent synthetic glioblastoma-like worlds with assumed signal quality `q` and assumed residual dependence `d`, is mean simulated utility higher for a controller that carries uncertainty over `d` than for an otherwise identical controller that fixes `d = 0`, and how does the difference change with the true `d` and `q`?"

- Negative: an interval at or below zero where `d` is nonzero and the signal is informative, or the `d = 0` controller wins.
- Artifact: a gain in the `q = 0.5` or null-channel controls.
- Inconclusive: a half-width above 14 days, an unadopted interval method, or a missing arm. On present evidence this is the expected category.
- Known weakness: the simulator generates data with `d`, so the `d = 0` controller is misspecified by construction. Any gain partly measures the cost of ignoring an assumption the simulator itself built in.

What the project can honestly claim today: one failed frozen prediction test (Track A), one exploratory imaging description on a spent cohort (Track B), a simulator with verified mechanics, a valid interval candidate, and a documented finding that the chosen precision target cannot be reached with available compute (Track C). It cannot claim patient drug resistance, clinical benefit, personalized treatment, or a validated integrated framework.

## 5. Verification record

- Scratch copy: `src/`, `tests/`, `configs/` copied outside the repository; 492 files compared byte-for-byte with the worktree after the last edit, 0 differences.
- Command (Python 3.14.4, `PYTHONPATH` set to the scratch `src`): `py -3.14 -m pytest tests/test_trackc_config_parity_extended.py tests/test_trackc_ge_betting.py tests/test_trackc_ge_betting_coverage_audit.py tests/test_trackc_integrated_sampler.py tests/test_trackc_world_summary.py tests/test_trackc_policy_benchmark.py tests/test_trackc_world_runner.py tests/test_track_b_transition_estimand.py -q -p no:cacheprovider`
- Result: 118 passed, 3 expected failures (the pinned drifts C9, C10, C11), 55.8 s.
- Screen command: `py -3.14 src/trackc_ge_betting_coverage_audit.py --processes 6`; 14.3 s; output copied unedited.
- Not run: the full test suite, any patient analysis, any policy outcome.
