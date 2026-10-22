# START HERE (for any AI agent or person). Updated 2026-10-21.

Read this file first. If another file disagrees with it, this file and `vault/PAPER_FINDINGS_LEDGER.md` win. This is an intentional repository instruction, not an empirical finding.

## 1. Current goal
Personal computational-biology research, not currently an ISEF/GARSEF submission. The user may clear or switch chats; keep this file current so the project can resume without chat history. User has 20 hours/week, mentors, a 16-core/34-GB desktop without a GPU, and Kaggle access for GPU runs. Check Kaggle quota, session limits, and dataset transfer needs before relying on a GPU; availability is not unlimited. Do not promise a prize or certainty of success.

## 2. Current project: "probabilistic RANO" (decision narrowed 2026-10-20)
Question: how do measurement noise and tool disagreement affect volumetric progression assessment and simulated PFS under RANO 2.0-inspired rules? RANO 2.0 uses context-specific reference scans: progression is assessed against the pretreatment baseline or the smallest measurement after treatment starts, with special rules in some settings; response is assessed against baseline. This project's simulations use a simplified running-minimum nadir. The 0.5 cm3 volume floor is our operational choice, not a RANO 2.0 rule. The project does not measure scan-rescan noise with its own pipeline.
- Plan: `paper/phd_level_plan_probabilistic_rano.md`. Problem and literature: `paper/problem_first_rano_noise.md`.
- Rules and declarations: `analysis_plan_v1.md` (Amendments 14 to 21 are the current ones). Results: `vault/PAPER_FINDINGS_LEDGER.md` (latest entry 2bx; next id 2by).
- Latest decision memo recommends NARROW, not a new theorem: the false-progression-by-visit-count calculation is routine and adds little beyond existing simulations. Any remaining contribution is applied sensitivity analysis, and its effect sizes depend on a noise value this project has not measured.
- Do not summarize the entire project as only an assumption audit: A15 is an exploratory real-data comparison of two tools' consecutive-scan classifications and noise proxies. It is not a direct scan-rescan noise measurement, and novelty is unverified.
- The project has tool disagreement, assumed noise in simulations, and a published scan-rescan anchor, but no scan-rescan result measured with its own pipeline. Abu Khalaf 2021 reports 19 of 20 participants overlap with RIDER; RIDER therefore is not independent calibration evidence. Treat 2–5 day change as potentially containing biology and systematic shift as well as noise.
- Final QIN feasibility memo (2026-10-20): `research_notes/final_feasibility_check_qin_path_2026-10-20_claude.md`. Decision: STOP the QIN path and do not write Amendment 22. None of the three reopening conditions is met: the Abu Khalaf supplement/within-subject SD scale is unavailable; QIN/dbGaP eligibility and cloud-GPU use are unverified and may require a qualified mentor; no segmentation method is validated on QIN's early-postoperative 5 mm scans. No access request or external message has been sent. The student may try the supplement in a normal browser and, if needed, send the draft author question in the memo; this is optional and does not reopen QIN by itself. Reopen only if measurement scale, legitimate access/compute, and a suitable segmentation route are all established. Any future analysis requires a declared and committed amendment before code runs, following rule 1.
- Latest memo: `research_notes/decision_memo_probabilistic_rano_applied_2026-10-20_claude.md` (may be ignored by git; key state is copied here). Preserve it and the other notes. Walovitch 2013 and whether Schmitt's slides cover nadir remain unresolved; do not cite/claim either without direct verification.
- Novelty is UNVERIFIED. Say "not found in our searches", never "first" or "revolutionary".

## 3. State of the work (checked against files on 2026-10-21)
| Item | State | Source |
|---|---|---|
| A15 test 1: consecutive-scan ratio class differs between the two tools | DONE. 54.5% of 121 eligible pairs (CI 45.6% to 63.8%), "material" | `output/rano_noise/rano_tool_dependence.json`, ledger 2br |
| A15 test 2: minimal detectable change vs the +40% line | DONE. Both proxy MDC95 estimates exceed line in smallest size band; estimates straddle line in two larger bands. These are sensitivity scenarios, not formal noise bounds. | same |
| A16 trial simulation (Aim 3) | DONE (simulation, assumed parameters). Noise makes PFS earlier and cuts power a lot (H-S1, H-S2 supported; H-S3, H-S4 NOT supported as declared) | `output/trial_sim/trial_sim_summary.json`, ledger 2bs |
| Calculator (`src/rano_calc.py`) and patient timeline viewer (`output/viewer_rano/viewer_rano.html`) | DONE (display and arithmetic only) | ledger 2bv, 2bw |
| Technical note | Draft audit reconciled 2026-10-21; external clinical/imaging review still needed before citation | `paper/technical_note_probabilistic_rano_a15_a21_draft_2026-10-20_claude.md` |
| Poster | Not a current goal | current scope decision |

Weakest point of A16: effect sizes depend on the assumed noise model; s=0.25 alone already gives -15 weeks. Weakest points of A15: the two tools are not a scan-rescan pair; 71.9% of pairs fell below the floor or had a missing volume; one centre; exploratory.

## 4. Rules that must not break
These are prospective workflow rules, not evidence that every historical run complied.
1. Declare a test in `analysis_plan_v1.md` and commit it BEFORE the code runs. Commit the code BEFORE the run. Run once. No re-runs with other thresholds.
2. Every number cites a file in `output/` or the ledger. No number from memory.
3. Report negative results unchanged. Label: simulation, assumed, post hoc, exploratory, interim.
4. Never hand-edit output JSON. Fix the code and re-run (only if the plan allows a re-run).
5. Paths come from `__file__`, not the working directory.
6. After a script runs, check the output file exists and its size is above zero.
7. Past work is NOT good just because it exists. Judge it on evidence.
8. H-2 and H-3 test sets are SPENT (Amendment 12). They cannot be locked tests again.
9. Commit format: `Script N: brief description`. Check `git status --short` first. Never `git add -A`. Push after each commit.

## 5. The older GBM work (Tracks A, B, C)
It is background, not the goal. Honest summary: `README_STATUS.md` (corrected after an audit). Important facts for anyone reading old files:
- Old files call work "first", "novel", "revolutionary", "validated". Those words are NOT supported. Ignore them.
- RL results: the 25% gain is on a simulated cohort; the real-test result (20/21) uses assumed kill rates. Both are simulation.
- Growth-model forecasts did not beat "no change" on real data (persistence baseline).
- Several drug and RL numbers use assumed parameters, not measured ones.
- Inventory checked 2026-10-19: 159 numbered top-level files in `src/` (maximum prefix 145, with collisions); 208 top-level and 217 recursive source `.py` files. `tests/` has 18 `test_*.py` files, 117 test functions, and conftest.py. File existence is not a safe-runnable or validation claim.

## 6. Which files to trust
| Status | Files |
|---|---|
| CURRENT | `START_HERE.md`, `CLAUDE.md`, `analysis_plan_v1.md`, `vault/PAPER_FINDINGS_LEDGER.md`, `paper/phd_level_plan_probabilistic_rano.md`, `paper/problem_first_rano_noise.md`, `paper/overview_for_review_rano.md` (dated review snapshot), `README.md` (new), `README_STATUS.md` (old-work status, corrected) |
| BACKGROUND (still true as records, wrong as goals) | `paper/new_project_options.md`, `paper/literature_search_digital_twin.md`, `vault/LIMITS_CURRENT.md`, `vault/claim_ledger.md`, `vault/audit_2026-10-18/*` |
| ARCHIVED (moved to `archive_old_gbm/`, same sub-paths; banner on each; do not use for goals or "first" claims) | all old GBM READMEs, summaries, plans, track results, `docs/*`, `output/*.md`, `vault/Script-*.md`, `vault/TRACK_*.md` |
| ARCHIVE (retained as records) | `vault/audit_2026-10-18/originals/*`, `.claude/skills/*`, `graphify-out/*` |

## 7. Next step
Aim 3 is done. Amendment 18 is DONE (script 142, ledger 2bt: bias appears from a change-SD of about 0.10; correlation shrinks it but does not remove it; flip rate stays 53% to 63% across floors and 54% against the nadir, n=68 pairs). Amendment 19 is DONE (script 143, ledger 2bu): with the change-SD taken from a published 2-day repeat-scan study (assumed conversion, large lesions, one study), the model still gives about 8 to 20 weeks earlier PFS and a large power loss. Calculator DONE (`src/rano_calc.py`, ledger 2bv; arithmetic only, not clinical). Timeline viewer DONE (`output/viewer_rano/viewer_rano.html`, ledger 2bw). Amendment 21 is DONE (script 145, ledger 2bx): A21's seven-value selected grid gives modeled single-scan median PFS shifts from -2.6 to -15.7 weeks at N=100 and 8-week scans. Across 11 distinct simulated values from A19 and A21, results also include -20.4 weeks at A19's s=0.342. The separately quoted -5.8 to -15.0-week range covers point-estimate conversion readings, whose largest shifts come from A19; it is not the full A21-grid range. Only 11 of 17 converted values were simulated. These are simulation results, not a measured project noise estimate.

Current phase: mentor gave positive overall feedback, per user; detailed comments and reviewer qualifications are not recorded. Desk-only dataset screen found QIN-GBM Treatment Response as strongest technical fit (double-baseline MRI, 54 patients, 33.5 GB), but its page requires a restricted license and TCIA's post-July-2025 controlled-data migration makes the current access route unclear. Kaggle use is unconfirmed. Open longitudinal alternatives do not provide clean scan-rescan noise estimates. No qualifying open repeat-scan tumor dataset has been verified in this preliminary search. See `paper/research_plan_next_phase.md` for sources and stop/go gates. No data download or new code until access, compute, and measurement-method requirements pass and a new amendment is committed. QIN remains stopped pending clarification; RIDER is restricted. Competition paperwork and poster are not current goals unless the user changes scope.

## 8. Known issues (from the Codex audit of 2026-10-19; see `paper/codex_audit_report.md`)
- Old scripts 15, 18, 20, 49 and 58 use working-directory-relative output paths (against the CLAUDE.md rule). Not fixed; they belong to the archived GBM work.
- `src/49_interactive_3d_dashboard.py` and `src/50_clinical_cdss_app.py` had merge-conflict markers (identical both sides); resolved 2026-10-19; static check on 2026-10-19 confirms no line-start conflict markers remain in either script. The earlier audit report describes the pre-fix state.
- Archived documents in `archive_old_gbm/` still contain unsupported claims ("first", "clinically viable", "validated"). Each has a banner. Do not quote them.
- The first audit sampled pre-2br ledger entries and old Markdown. Round 2 completed its full claim checklist (`paper/codex_audit2/SUMMARY.md`); corrections remain in progress in `paper/codex_audit2/fix_progress.md`. Audit completion is not scientific validation.
- Scripts 140 to 145 were not re-run by the auditor; outputs were compared statically.
- `output/drug_gating_report.md` and similar old reports were moved to `archive_old_gbm/output/`; scripts that name them will not find them.
