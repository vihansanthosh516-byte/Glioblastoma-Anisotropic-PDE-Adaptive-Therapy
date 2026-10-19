# Audit notes: documents, outputs, and repository hygiene

Read-only audit of the current RANO-facing documents, cited RANO JSON/CSV outputs, root and `paper/` inventories, and tracked file sizes. Older archived claims and old-work output verification are assigned to the old-work audit. Findings below use repository-relative file and line numbers.

## F. Claims in current documents

| Finding | File:line | Evidence and check |
|---|---|---|
| Stale completion status says calculator and viewer are not started. | `START_HERE.md:20` | Reads “Calculator, patient timeline viewer, paper, poster \| NOT started.” The ledger entries 2bv and 2bw say calculator and viewer are done, and `README.md:19-21` links their code/output. |
| Stale next-step overview says calculator and timeline viewer are not done. | `paper/overview_for_review_rano.md:36-37` | Reads “Identifiability proof, calculator, patient timeline viewer, paper, poster, clinician review.” Calculator and timeline viewer are complete per ledger 2bv/2bw and `README.md:19-21`. Identifiability proof, paper, poster and clinician review remain open. |
| Stale plan status says A4 calculator and viewer are not done. | `paper/phd_level_plan_probabilistic_rano.md:3` | Reads “A1 proof, A4 calculator and viewer are not done.” A4 tools are done (ledger 2bv/2bw; `README.md:19-21`). The A1 proof still appears incomplete. |
| Stale research brief claims nothing is declared or run, and describes A15 analyses as future. | `paper/problem_first_rano_noise.md:3,27-28` | Line 3 reads “Nothing here is declared or run.” Lines 27-28 describe a planned viewer and tests; A15 tests and timeline viewer already exist (ledger 2br/2bw; `README.md:19-21`). Treat this file as a historical draft or update its status. |
| A21 summary omits the largest reported power loss among its evaluated point estimates. | `paper/overview_for_review_rano.md:29` | It reports power change “(-0.07 to -0.46).” `output/conversion_audit/conversion_audit_summary.json` has `N100_D8_single_scan["0.262"].E2 = -0.703` and `E2_hi = -0.6827`; the ledger 2bx text acknowledges “-0.07 to -0.46 or more.” Wording understates the saved output. |
| 0.5 cm3 floor is explicitly described as a project choice in current summary materials; no violation found in these files. | `README.md:6`, `START_HERE.md:9`, `paper/overview_for_review_rano.md:6,23`, `paper/phd_level_plan_probabilistic_rano.md:6`, `paper/problem_first_rano_noise.md:6` | Each says the floor is operational/ours and not part of RANO 2.0. The cited `output/rano_noise/rano_tool_dependence.json` records a 500 mm3 exclusion threshold. |

No current RANO-facing doc sentence claiming “first,” “novel,” “validated,” or “proves” was found as an affirmative project result; the current docs consistently mark novelty unverified and the calculator/viewer non-clinical. The old-document sweep is delegated separately.

## I. Repo hygiene and doc conflicts

| Finding | File:line | Evidence and check |
|---|---|---|
| `START_HERE.md` inventory classifies only selected current/background/archived Markdown; several tracked `paper/` notes have no explicit status there. | `START_HERE.md:43-49` | `git ls-files 'paper/*.md'` shows `paper/overview_for_review_rano.md`, `chatgpt_review_response.md`, `codex_audit_prompt.md`, `research_prompt_v2.md`, and `research_v2_quick_pass.md` beyond the current/background list. Some are plainly historical/review artifacts, but status is ambiguous. |
| The research brief and independent-review overview conflict with the current work state. | `paper/problem_first_rano_noise.md:3,27-28`; `paper/overview_for_review_rano.md:36-37`; `START_HERE.md:20,52` | Both claim completed calculator/viewer tasks remain undone; see F. |
| Multiple top-level and `paper/` work products look historical or duplicated, but are not all labeled in `START_HERE.md`. | `START_HERE.md:43-49` | Tracked root docs include `README_STATUS.md` and `CLAUDE.md`; `paper/` contains two research prompt/quick-pass files and older project-option/literature notes. No deletion recommended without owner review. |

## M. Cross-document or source-number discrepancies

| Finding | File:line | Evidence and check |
|---|---|---|
| A21 power-loss range differs across summary and ledger, and the summary omits the largest evaluated value. | `paper/overview_for_review_rano.md:29`; `vault/PAPER_FINDINGS_LEDGER.md:130` | Overview says -0.07 to -0.46. Ledger says -0.07 to -0.46 “or more.” Saved JSON contains `E2 = -0.703` at s=0.262. Cite JSON and report the full evaluated range (-0.029 to -0.703 for the listed point-estimate grid) or describe the bound precisely. |
| A15’s 54.5% estimate has two slightly different upper CI limits in separate saved analyses. These are not necessarily a contradiction because they are separate bootstrap runs. | `README.md:11`; `output/rano_noise/rano_tool_dependence.json`; `output/rano_sensitivity/rano_sensitivity_summary.json` | A15 CI upper endpoint is 0.638095 (63.8%). A18 floor-500 endpoint is 0.635419 (63.5%), same 121 pairs and rate. Summaries should identify which bootstrap/result they cite, rather than imply one shared CI. |

Other current README point estimates inspected matched the cited JSON at displayed precision: A15 rate/pairs (0.54545, 121); A18 nadir rate/sample (0.54412, 68 pairs/20 patients); A19 PFS range (-7.644 to -20.44 weeks; power changes -0.1995 to -0.82); A21 PFS point-estimate range (-5.768 to -14.988 weeks); A16 power values (0.969, 0.288, 0.054). Files opened: `output/rano_noise/rano_tool_dependence.json`, `output/rano_sensitivity/rano_sensitivity_summary.json`, `output/anchored_sim/anchored_sim_summary.json`, `output/conversion_audit/conversion_audit_summary.json`, `output/trial_sim/trial_sim_summary.json`, `output/viewer_rano/viewer_rano_info.json`; each exists and has nonzero size. No cited RANO JSON/CSV path checked here was missing.

## N. Dead files and large tracked files

| Finding | File:line | Evidence and check |
|---|---|---|
| One tracked source data file exceeds 10 MB. | `data/tcga_gbm_expression.tsv` | `git ls-files` plus filesystem size check: 24.42 MB. It is a large tracked data artifact. |
| Empty run logs and a session transcript exist at repository root; likely abandoned run artifacts. | `run_hash_0.log`, `run1.log`, `terminal_session.txt` | Root inventory: both logs are 0 bytes; transcript is 970 bytes. `git ls-files '*log*' '*err*'` did not list these root files, so they are untracked/local. |
| Multiple untracked cache/log/run directories appear in current Git status. | `git status --short` (no line number; worktree status) | Includes `kaggle_run/*`, `output/*/cache*`, `output/pde_manifest/*.log.err`, `output/shape_skill/`, and a large viewer artifact. These are untracked at audit time; distinguish generated caches from deliverables before cleanup. |
| A small untracked one-off Python patch helper sits in the paper folder. | `paper/_apply7.py` | Present in `git status --short` as untracked; no references found in tracked docs during targeted search. Looks like a temporary edit helper. |

Tracked files above 10 MB: `data/tcga_gbm_expression.tsv` (24.42 MB). No other tracked file exceeded the threshold in this checkout.

## Output reference sweep

| Finding | File:line | Evidence and check |
|---|---|---|
| The tracked-Markdown sweep found 17 concrete missing legacy output references, all inside an archived planning document; one additional `output/file.json` is an example path in the run instructions. | `archive_old_gbm/paper/GRAND_PLAN.md:90-96,102-103,113-119,125-126,133-135`; `CLAUDE.md:104` | Scanned all 104 tracked Markdown files and 134 unique output JSON/CSV references against the checkout and `archive_old_gbm/`. Missing refs: `output/pde_manifest/results_final.json`, `output/predict_gbm/results.json`, `output/forecastability_measured.json`, `output/pde_manifest/first_vs_later.json`, `output/volume_forecast/results.json`, `output/selective_forecast.json`, `output/power_sim.json`, `output/ipw_sensitivity.json`, `output/pde_manifest/pilot_1mm.json`, `output/solver_verification_v2.json`, `output/identifiability.json`, `output/patient_calibrated.json`, `output/kill_rate_fit.json`, `output/rho_continuous.json`, `output/trackA_gbmap.json`, `output/trackA_core_vs_periphery.json`, `output/trackA_balanced.json`. `CLAUDE.md:104` uses `output/file.json` as an illustrative example, not a claim that it exists. All current RANO output references inspected exist. |

## Limits

This note audits current RANO-facing Markdown and performs a tracked-Markdown output-reference sweep, but does not audit every archived/old-work claim. Those claims belong to the separate old-work audit. “No missing cited RANO output found” applies to the current summaries and ledger entries 2br-2bx inspected here.
