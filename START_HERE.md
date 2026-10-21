# START HERE (for any AI agent or person). Updated 2026-10-19.

Read this file first. If another file disagrees with it, this file and `vault/PAPER_FINDINGS_LEDGER.md` win. This is an intentional repository instruction, not an empirical finding.

## 1. Current goal
One author. Greater Austin Regional Science and Engineering Fair (GARSEF; tentative deadline February 2027; official deadline, dates, categories and forms NOT verified). Goal: a new, strong computational-biology project with real mathematics, backed by solid science. A win is not predictable. Never promise one.

## 2. Current project: "probabilistic RANO"
Question: are the RANO 2.0 volume lines (progression at +40%, partial response at -65%; progression is judged against the NADIR, response against baseline). The 0.5 cm3 volume floor used in our analysis is OUR operational choice. It is NOT in the RANO 2.0 text (checked 2026-10-18 in PMC11630866) above the measurement noise of automated segmentation in glioblastoma, and what does the noise do to trial endpoints (progression-free survival, PFS)?
- Plan: `paper/phd_level_plan_probabilistic_rano.md`. Problem and literature: `paper/problem_first_rano_noise.md`.
- Rules and declarations: `analysis_plan_v1.md` (Amendments 14 to 21 are the current ones). Results: `vault/PAPER_FINDINGS_LEDGER.md` (latest entry 2bx; next id 2by).
- Novelty is UNVERIFIED (limited snippet-level searches; no retained count log). Say "not found in our searches", never "first" or "revolutionary".

## 3. State of the work (checked against files on 2026-10-18)
| Item | State | Source |
|---|---|---|
| A15 test 1: consecutive-scan ratio class differs between the two tools | DONE. 54.5% of 121 eligible pairs (CI 45.6% to 63.8%), "material" | `output/rano_noise/rano_tool_dependence.json`, ledger 2br |
| A15 test 2: minimal detectable change vs the +40% line | DONE. Both proxy MDC95 estimates exceed line in smallest size band; estimates straddle line in two larger bands. These are sensitivity scenarios, not formal noise bounds. | same |
| A16 trial simulation (Aim 3) | DONE (simulation, assumed parameters). Noise makes PFS earlier and cuts power a lot (H-S1, H-S2 supported; H-S3, H-S4 NOT supported as declared) | `output/trial_sim/trial_sim_summary.json`, ledger 2bs |
| Calculator (`src/rano_calc.py`) and patient timeline viewer (`output/viewer_rano/viewer_rano.html`) | DONE (display and arithmetic only) | ledger 2bv, 2bw |
| Paper, poster | Author-reported not started as of 2026-10-18 | plan; no draft-status audit |

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
Aim 3 is done. Amendment 18 is DONE (script 142, ledger 2bt: bias appears from a change-SD of about 0.10; correlation shrinks it but does not remove it; flip rate stays 53% to 63% across floors and 54% against the nadir, n=68 pairs). Amendment 19 is DONE (script 143, ledger 2bu): with the change-SD taken from a published 2-day repeat-scan study (assumed conversion, large lesions, one study), the model still gives about 8 to 20 weeks earlier PFS and a large power loss. Calculator DONE (`src/rano_calc.py`, 7 of 7 declared tests pass, ledger 2bv; arithmetic only, not clinical). Timeline viewer DONE (`output/viewer_rano/viewer_rano.html`, ledger 2bw). Amendment 21 is DONE (script 145, ledger 2bx): the published repeatability can be converted to a noise level in three ways (change-SD 0.079 to 0.342); the PFS bias direction holds on the seven selected simulated values (not every one of the 17 converted values), with modeled single-scan median PFS bias about 5.8 to 15.0 weeks at N=100 and 8-week scans, and at the lowest reading the two-scan rule keeps the bias under a month. Next-step plan recorded 2026-10-19: paper draft and poster. Check GARSEF 2027 handbook (dates, categories, forms) and ask a clinician or faculty member to review the framing.

## 8. Known issues (from the Codex audit of 2026-10-19; see `paper/codex_audit_report.md`)
- Old scripts 15, 18, 20, 49 and 58 use working-directory-relative output paths (against the CLAUDE.md rule). Not fixed; they belong to the archived GBM work.
- `src/49_interactive_3d_dashboard.py` and `src/50_clinical_cdss_app.py` had merge-conflict markers (identical both sides); resolved 2026-10-19; static check on 2026-10-19 confirms no line-start conflict markers remain in either script. The earlier audit report describes the pre-fix state.
- Archived documents in `archive_old_gbm/` still contain unsupported claims ("first", "clinically viable", "validated"). Each has a banner. Do not quote them.
- The first audit sampled pre-2br ledger entries and old Markdown. Round 2 completed its full claim checklist (`paper/codex_audit2/SUMMARY.md`); corrections remain in progress in `paper/codex_audit2/fix_progress.md`. Audit completion is not scientific validation.
- Scripts 140 to 145 were not re-run by the auditor; outputs were compared statically.
- `output/drug_gating_report.md` and similar old reports were moved to `archive_old_gbm/output/`; scripts that name them will not find them.
