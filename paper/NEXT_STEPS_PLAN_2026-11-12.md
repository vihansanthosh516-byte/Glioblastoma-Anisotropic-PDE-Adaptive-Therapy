# Complete next-steps plan (written 2026-11-12 for handoff to Claude or Codex)

Read this file after `START_HERE.md`. It records the current state and next safe work. On resume, first run `git log --oneline -8` and `git status --short`, and check that the files named here exist.

> **Track B status correction (2026-11-12):** ACRIN 6677 area forecasting was completed and failed its frozen rule. Its separate sRCBV linkage audit found 15 eligible patients against a frozen minimum of 25 and stopped before reading measurement magnitudes. Sections 2–8 below preserve the old ACRIN procedure for audit history only; do not execute them. See `paper/trackB_acrin6677_result_2026-11-11.md` and `paper/trackB_acrin6677_advanced_mri_linkage_result_2026-11-11.md`.

## 0. Where things stand (verified on disk 2026-11-12)

| Part | State | Record |
|---|---|---|
| Track A A5 | Failed frozen rule (4 of 6 donors, rule 5 of 6). Stays failed. | `output/trackA_yu2020/result.json` |
| Track B Duke | Exploratory, cohort spent. | `paper/gap_status_2026-11-09.md` B1 |
| Track B LUMIERE | No-go (7 of 28 passed QC). | `paper/trackB_lumiere_registration_nogo_2026-11-11.md` |
| **Track B ACRIN 6677** | Area forecast completed and failed. sRCBV exact-link feasibility audit stopped at 15/25 before measurement magnitudes. No further ACRIN model is active. | `paper/trackB_acrin6677_result_2026-11-11.md`, `paper/trackB_acrin6677_advanced_mri_linkage_result_2026-11-11.md` |
| Track C | Confirmatory evaluation **stopped**; 14-day target unchanged (user decision 2026-11-10). Do not reopen or suggest changing it. | `START_HERE.md` |
| Integration | No patient-level link among biology, MRI, and measured drug response exists in open data. | `paper/data_request_plan_2026-11-12.md` |
| HGCC bridge | Withdrawn: no measured TMZ response found. | same file |
| Data requests | User allowed them, then stated they cannot ask anyone. Drafts exist, unsent, parked. | `paper/data_request_drafts_2026-11-12.md` |

Data on disk: `data/acrin6677/ACRIN-DSC-MR-Brain-Full-Clinical.zip` (SHA-256 `584a17f37017e5a00962a1677cc4640a6073f79ec1ae905f639e0867d9866157`), extracted to `data/acrin6677/x/`. Folder is git-ignored. Licence CC BY 4.0 (TCIA record).

## 1. Hard rules (do not break)

1. Never tune, relax, or re-choose a test after seeing results. The user asked once to "bend the test"; this was declined and stays declined (project rule: never fabricate or tune toward a positive result).
2. The original pre-analysis ACRIN restriction is historical. The sRCBV feasibility rule now blocks further ACRIN measurement reads and model fitting on this cohort. Do not relax its 25-patient minimum.
3. The old ACRIN instructions in sections 2–8 are historical. No ACRIN 6677 measurement or model analysis is active.
4. Track C stays stopped. No policy outcomes.
5. Report whatever the result is: success, failure, or inconclusive. Put the weakest point next to the result.
6. Run Python that reads downloaded files with `python -I`. Scripts live in `src/`, never inside `data/acrin6677/`.
7. Paths in scripts derive from `__file__` (see `CLAUDE.md`).
8. Commit and push after each verified step. Check `git status --short` first. Never `git add -A`.
9. Subagents use `model: "haiku"`.

## 2. Step 1: get the user's answer on the S3 amendment (pending)

The question open in chat: "Write a dated amendment that takes the volume unit from the trial paper?" If the user says yes (or already did), go to Step 2. If the user says stop, write `paper/trackB_acrin6677_nogo_<date>.md` and jump to Section 9.

## 3. Step 2: S3 amendment (documentation only, before any value)

1. Find the unit of the central-reader 3D enhancing volume in the trial's own publications. Try in order:
   - Boxerman et al. 2013, Neuro Oncol 15:945-954, doi 10.1093/neuonc/not049 (PMC full text).
   - Boxerman et al. 2018/2019 ACRIN 6677 volumetric papers; Schmainda et al. 2015 (Neuro Oncol 17:1148) ACRIN 6677 rCBV paper.
   - The TCIA collection page "Detailed Description".
2. Quote the sentence (under 15 words) and the URL.
3. Write `paper/trackB_acrin6677_amendment1_<date>.md` containing:
   - the quote and source;
   - the unit (expected cm^3 or mm^3; do not assume);
   - the constant `c` for `y = log(v + c)`: **one tenth of one unit of the documented scale if the paper states no measurement resolution** (fix this rule now, here, so it is not chosen later: `c = 0.1` in the documented unit);
   - exact variable names: patient `cn`; time `S0e2d` (days); reader `s0e4deid`; measurement `s0e33`; visit label `s0e1`;
   - visit rule: a "measured visit" is a distinct `S0e2d` day with at least one non-missing numeric `s0e33`; reader values on the same day are averaged; visits ordered by day; first three measured visits form the window; the patient is eligible only if day1 < day2 < day3;
   - rows with `s0e32` = "No" (3D: no lesions present) and a missing volume: treated as volume 0 **only if** the amendment states it; decide from the dictionary meaning, not from data. Recommended and fixed here: lesion-absent with missing volume counts as 0; uninterpretable image quality (`s0e24` = 1) counts as missing;
   - a statement that no measurement value was read before the amendment commit.
4. If no publication states the unit: write the no-go record instead. Do not guess.
5. Commit and push. This commit is the amendment freeze.

## 4. Step 3: code and synthetic tests (still no patient values)

Create:

- `src/track_b_acrin6677_forecast.py`: pure functions only.
  - `predict_b1(y1,y2)`, `predict_b2(t1,y1,t2,y2,t3)`, `predict_b3(...)` (damping 0.5), `predict_m(...)` with `tau = 56` days.
  - Model M: `dy/dt = g - k*exp(-t/tau)`, t measured from visit 1. Integrated: `y(t) = y1 + g*(t - t1') - k*tau*(1 - exp(-(t)/tau))` with `t` = days since visit 1. Two history points give one increment (visit 1 to visit 2), so two parameters are **not** identified by two points. The protocol says "solved exactly from the two history increments"; with only two history visits there is one increment. **This is a protocol defect found at planning time, before any value was read.** Resolve it in the Step 2 amendment, not later, with this fixed rule: set `g = 0` when `y2 < y1` and fit `k` from the single increment (pure decaying treatment effect, then plateau); set `k = 0` when `y2 >= y1` and fit `g` (pure growth, equals B2). State plainly that M then differs from B2 only for patients whose volume fell between visits 1 and 2, where M predicts a flattening decline and B2 predicts continued decline.
  - `score(pred, obs) = abs(pred - obs)`.
- `src/track_b_acrin6677_run.py`: loads S0 from both folders, builds windows by the amendment rules, writes `output/track_b_acrin6677_<date>/result.json` and `per_patient.csv`. Has a `--dry-run` flag that prints only the eligible count.
- `tests/test_track_b_acrin6677_forecast.py` (protocol section 6): (a) M recovers parameters on noise-free M curves in both branches; (b) on pure log-linear rising data M equals B2; (c) finite outputs for flat, rising, falling histories; (d) scoring matches a hand-built three-patient table; (e) window builder on a synthetic S0-like table: reader averaging, day ordering, lesion-absent rule, exclusion when fewer than three visits.

Run tests in a scratch copy (see `CLAUDE.md` Testing; set `PYTHONPATH` to scratch `src`). All must pass. Commit code and tests with the test output pasted in `paper/trackB_acrin6677_test_record_<date>.md`. Record the code SHA-256.

## 5. Step 4: the one patient run

1. Confirm `git status` is clean for the three source files and that their hashes match the test record.
2. `python -I src/track_b_acrin6677_run.py --dry-run` : eligible count only. If fewer than 25, stop and write a no-go.
3. `python -I src/track_b_acrin6677_run.py` : once. No reruns with changed settings.
4. Verify on disk (`rules/verify.md`): `result.json` exists, size > 0, fields present: n, per-method median error with bootstrap 95% CI (10,000 resamples, seed 20261112), best baseline name, Wilcoxon p, median paired difference, win count, verdict.

Decision rule (frozen, protocol section 5): success = median paired difference favours M **and** two-sided exact Wilcoxon p < 0.05 against the best of B1/B2/B3; failure = a baseline is significantly better or the median favours it; otherwise inconclusive.

Honest expectation: M equals B2 for patients whose volume rose, so the test is decided only by patients whose volume fell after bevacizumab. The result can go either way.

## 6. Step 5: report

Write `paper/trackB_acrin6677_result_<date>.md` in the project's 4-question format (memory: result-report-format): what was tested, what was found, weakest point, what next. Include all three baselines, the amendment, the planning-time protocol defect, pseudo-response, single window, and that novelty is unchecked. Update `START_HERE.md`, `tracks/track-b/README.md`, and `vault/PAPER_FINDINGS_LEDGER.md`. Commit and push.

- If **failure or inconclusive**: per project rule, launch one haiku improvement subagent to propose the next pre-registered idea on **new** data or a clearly labelled exploratory analysis. Do not re-test this cohort as confirmatory.
- If **success**: run the prior-art search (Section 7) before any novelty wording.

## 7. Step 6: prior-art search for the ACRIN result

Search primary literature for: forecasting tumor volume after bevacizumab in recurrent GBM; ACRIN 6677 / RTOG 0625 secondary analyses; growth-kinetic models under anti-angiogenic therapy (for example Swanson-group proliferation-invasion models; Ellingson volumetric response work; tumor growth inhibition models in oncology pharmacometrics, for example Claret-type TGI models). Write `paper/trackB_acrin6677_prior_art_<date>.md` saying what exists, how it was checked, and what is unverified. Expect that tumor-growth-inhibition models are well established: the honest novelty claim is narrow (this cohort, this pre-registered comparison).

## 8. Exploratory extras on ACRIN 6677 (only after Step 5, each labelled exploratory)

These use the same open file and need a short dated note before each one:

1. Same comparison on 2D area (`s0e23`) and T2 volume (`s0e43`).
2. Later windows (visits 2,3 to 4).
3. By treatment arm (arm variable: find in form A0 dictionary; not yet identified).
4. Early volume change versus survival (death days in form DS, `DSe5d`). Note: Boxerman 2013 already published early progression versus survival on this trial, so this is confirmatory of known work, not new.
5. **Completed feasibility check:** the SC sRCBV exact-link audit found 15 eligible patients against the frozen 25-patient minimum. It stopped before measurement magnitudes; no advanced-MRI model was run. See `paper/trackB_acrin6677_advanced_mri_linkage_result_2026-11-11.md`. The old suggestion to substitute DE or MS is withdrawn.

## 9. Integration status and next route

No patient-level link among Track A expression, Track B MRI, and measured drug response exists in the current evidence. The ACRIN 6677 area forecast failed its frozen rule; its sRCBV analysis stopped at 15 of 25 before reading values. It does not validate the Track C response curve. Track C policy evaluation remains stopped under the user's unchanged 14-day target.

- **Track A to drug response:** Ntafoulis 2023 has RNA-seq for at most 24 cultures with measured in-vitro TMZ response. This could support a small descriptive association only; it cannot establish patient drug resistance or repair A5's failed rule.
- **Track A to Track B:** Ivy GAP has region RNA for at most 10 tumors and images require controlled access. It cannot meet the 25-patient gate.
- **Integrated claim:** any future A/B-to-C connection remains a simulation unless a matched patient cohort supplies the required measurements. Do not describe simulated treatment benefit as clinical benefit.

## 10. Fresh-cohort screening after ACRIN no-go

Keep all cohort checks metadata-only until an eligible source, access terms, and exact repeated-measure linkage are established. Current candidate screen:

1. **ACRIN 6686 / RTOG 0825: strongest request-gated lead.** The trial protocol scheduled serial MRIs and central review. A published central-reader analysis used 284 patients with two timepoints only; it does not establish how many have three or more linked quantitative visits. The advanced-MRI subset enrolled 44 patients, with smaller per-timepoint evaluable counts. Requesting data requires the ECOG-ACRIN Data Access Committee; access, output terms, complete visit counts, and overlap with prior project cohorts remain unresolved. This is a candidate, not an eligible cohort. Sources: `https://www.acr.org/Clinical-Resources/Publications-and-Research/Clinical-Research/acrin-legacy-trials`, `https://pmc.ncbi.nlm.nih.gov/articles/PMC6120359/`, `https://pmc.ncbi.nlm.nih.gov/articles/PMC7906067/`.
2. **BraTS Longitudinal:** a published benchmark reports 164 glioma patients with three to five scans each, but this project has prior BraTS exposure. Exact patient overlap, disease-specific eligible count, and dataset-use terms must be checked before any fresh-confirmation claim. Treat any reused or exposed subjects as exploratory. Source: `https://www.synapse.org/#!Synapse:syn28740011/wiki/70861`.
3. **Burdenko-GBM-Progression:** 180 patients and serial imaging are reported, but images require controlled access, follow-up GTV labels cover only an unknown subset, and this project already used its clinical CSV. It is not fresh confirmation. Source: `https://www.cancerimagingarchive.net/collection/burdenko-gbm-progression/`.
4. QIN-GBM Treatment Response is a repeatability study with double-baseline scans; Brain-Tumor-Progression has 20 patients and two exams; QIN-BRAIN-DSC-MRI has 49 subjects and 52 studies. These do not meet the current three-visit, 25-patient forecast gates.

The next research task is to finish a source-checked, request-ready feasibility screen for ACRIN 6686 and BraTS Longitudinal. Do not download, contact anyone, read patient rows, or adopt a protocol until metadata availability, access terms, and overlap gates are clear. Do not lower the cohort threshold or change outcomes to obtain a positive result.

## 11. External actions and records

- ACRIN legacy data requests must go through the ECOG-ACRIN Data Access Committee. No request has been sent. A qualified sponsor may be required; the public access page does not establish that a trainee can apply independently.
- The Duke result remains exploratory and its accepted license/share terms are not saved. Do not share Duke-derived summaries or weights until terms are verified.
- The Zenodo access token was pasted into chat earlier. Revoke it and create a replacement privately; never store it in the repository or chat.
- Track C's 14-day target and stop decision remain unchanged. Open code-setting decisions do not authorize policy outcome runs.

## 12. Claim boundary

No plan can guarantee a positive result or uniqueness. New analysis must use a predeclared question, cohort, endpoint, comparator, split, uncertainty method, and failure rule. Report failed or inconclusive results as observed. A better result cannot be created by relaxing eligibility, changing outcomes after inspection, or presenting simulation as patient evidence.
