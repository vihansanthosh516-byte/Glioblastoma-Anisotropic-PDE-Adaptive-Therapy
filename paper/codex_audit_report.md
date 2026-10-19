# Independent repository audit

## Summary

| Area | Wrong | Stale | Overclaim | Highest-risk finding |
|---|---:|---:|---:|---|
| A-D, G: current RANO | 0 headline mismatches | 1 | 1 | A21 calls seven simulated values a union, but conversion table has 18 distinct values. |
| E, J: tests and environment | 0 test failures in final scratch run | 0 | 0 | `108 passed`; two optimizer warnings and one convergence warning. |
| F, I, M, N: docs and hygiene | 1 numeric summary range | 4 stale status statements | 0 affirmative current overclaims found | A21 overview omits saved power loss `-0.703`. |
| H: leakage and reuse | 0 confirmed pre-freeze reads | 0 | 0 | H-2/H-3 sets are explicitly spent and reused as development after Amendment 12. |
| K-L: older claims and code | 2 syntax failures | Several | Many historical claims | Legacy claims remain in archived docs; old-code audit is sampled, not exhaustive. |

Audit notes: [RANO](</C:/Users/vihan/20206 science fair/paper/codex_audit/rano.md>), [old tracks](</C:/Users/vihan/20206 science fair/paper/codex_audit/old_tracks.md>), [docs and outputs](</C:/Users/vihan/20206 science fair/paper/codex_audit/docs_outputs.md>).

## A. Number check

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| Ledger 2br–2bx headline numbers match cited outputs to reported precision. A15 pair CSV independently gives 121 rows and flip rate 0.54545; output CI is 0.4559–0.6381. | `vault/PAPER_FINDINGS_LEDGER.md:124-130`; `output/rano_noise/rano_tool_dependence.json` | Opened summaries and cited CSVs; compared values at displayed precision. | Low |
| No missing cited output found for entries 2br–2bx. | `vault/PAPER_FINDINGS_LEDGER.md:124-130` | Checked JSON/CSV summary paths in the entries; all opened with nonzero content. | Low |

## B. Order check

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| Amendment declaration → code → result order passes for Amendments 15–21. | Git history; examples `51fd053 → ad9c55c → 945e2ab`, `2d3e6a4 → 48c5933 → 131495f` | Compared commit subjects and chronology with `git log`. | Low |

## C. Declaration vs code

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| A21 simulation grid does not cover stated conversion-value union. Script runs seven rounded `s` values; conversion table has 18 distinct values, including 0.291 and 0.342. PU1/PU2 apply only to the seven-point grid checked. | `analysis_plan_v1.md:624`; `src/145_conversion_audit.py:43`; `output/conversion_audit/conversion_table.csv` | Compared code grid with all conversion rows and independently recalculated R1/R2/R3 arithmetic. | High |
| Other inspected formulas, seeds, estimands, thresholds, calculator rules, and pass rules match their amendments. | `src/140_rano_tool_dependence.py:57-69`; `src/141_trial_sim.py:42-115`; `src/142_rano_sensitivity.py:94-101`; `src/145_conversion_audit.py:47-53`; `src/rano_calc.py:28-78` | Static comparison with Amendments 15–21. No scripts 140–145 were run. | Medium |

## D. Prediction check

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| Saved summaries support H-S1/H-S2, reject H-S3/H-S4, and support PA1/PA2/PB1/PD1, PX1–PX4, PU1–PU4 as ledger states them. | `output/trial_sim/trial_sim_summary.json`; `output/rano_sensitivity/rano_sensitivity_summary.json`; `output/anchored_sim/anchored_sim_summary.json`; `output/conversion_audit/conversion_audit_summary.json` | Recomputed support status from saved booleans and declared CI/materiality criteria. PU results are limited by grid issue in C. | Medium |

## E. Tests

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| Requested calculator tests pass: 7 passed. | `tests/test_rano_calc.py` | Ran `python -m pytest tests/test_rano_calc.py -q`; output: `7 passed in 0.68s`. | Low |
| Full suite passes in Git-backed scratch clone: 108 passed, 13 warnings. | `tests/` | Ran `python -m pytest tests/ -q` in `C:\Users\vihan\AppData\Local\Temp\science-fair-audit-clone-20261019`; warnings include unknown `disp` solver option and one sklearn convergence warning. | Low |

## F. Claims in current docs

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| Status incorrectly says calculator and patient timeline viewer are not started. | `START_HERE.md:20` | Ledger 2bv/2bw and `README.md:19-21` show both complete. | Medium |
| Overview and plan still list calculator/viewer as unfinished; brief says nothing is declared or run and describes completed work as future. | `paper/overview_for_review_rano.md:36-37`; `paper/phd_level_plan_probabilistic_rano.md:3`; `paper/problem_first_rano_noise.md:3,27-28` | Compared claims with ledger 2br/2bv/2bw and current README. | Medium |
| A21 overview gives power range `-0.07 to -0.46`, omitting saved `E2=-0.703` at `s=0.262`. | `paper/overview_for_review_rano.md:29`; `output/conversion_audit/conversion_audit_summary.json` | Compared overview range with the JSON point-estimate grid and ledger 2bx. | High |
| No current RANO doc calls 0.5 cm³ a RANO 2.0 rule. It is described as the project’s operational choice. | `README.md:6`; `START_HERE.md:9`; plan/overview/brief | Searched and read current RANO-facing docs; checked 500 mm³ output threshold. | Low |

## G. Science check

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| KM scan-grid risk sets and log-rank score/variance are internally consistent for discrete scan times. PFS uses prior-scan nadir and a simplified two-scan confirmation against the same nadir. | `src/141_trial_sim.py:42-60,63-76,94-115` | Traced event and censoring code against Amendment 16. This is a simulation rule, not full RANO adjudication. | Medium |
| AR(1) recurrence uses stationary marginal SD; patient-cluster bootstrap resamples patients; conversion arithmetic matches declared equations. | `src/142_rano_sensitivity.py:94-101`; `src/140_rano_tool_dependence.py:57-69`; `src/145_conversion_audit.py:47-53` | Static code inspection and independent arithmetic. | Low |
| A15 consecutive-scan ratio is called a RANO 2.0 volumetric class in its original ledger entry, though PD is nadir-referenced and PR baseline-referenced. A18 later acknowledges this. | `src/140_rano_tool_dependence.py:10`; `vault/PAPER_FINDINGS_LEDGER.md:124,126`; `analysis_plan_v1.md:550` | Compared A15 implementation with the A18 correction. | Medium |
| Simulation omits non-volume progression criteria, human segmentation approval, treatment-era rules, and full confirmation workflow. | `src/141_trial_sim.py:42-60,94-115`; `analysis_plan_v1.md:509-516,578-585,624-630` | Compared modeled event process and declared limitations. | Medium |

## H. Leakage and reuse

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| H-2 loader blocks test recurrence access until both freeze config and `h2-frozen` tag exist. H-3 scripts require freeze config/tag and check frozen input hashes. | `src/predictgbm_io.py:81-90`; `analysis_plan_v1.md:324,337`; `configs/h3_frozen.yaml:1` | Inspected guard code, configuration, tags, and declaration/result chronology. | Medium |
| H-2 and H-3 sets are no longer locked tests. Amendment 12 explicitly relabels LUMIERE/RHUH and UCSF as post hoc development data. | `analysis_plan_v1.md:402-408,450` | Read reuse declaration and compared with later analysis declarations. | High |
| No pre-declaration outcome read was confirmed in inspected source/history. Direct external/manual access before freezes cannot be established from repository files. | `analysis_plan_v1.md:318-337`; commits `cd20f0d`, `fb1a270` | Checked stated freeze chronology and source guards; this does not prove external access history. | Medium |

## I. Repo hygiene

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| Several tracked `paper/` review and research notes lack clear current/background/archive status in `START_HERE.md`. | `START_HERE.md:43-49` | Compared inventory with tracked Markdown files. | Low |
| 17 missing legacy JSON/CSV references occur in archived `GRAND_PLAN.md`; current RANO cited outputs checked here exist. | `archive_old_gbm/paper/GRAND_PLAN.md:90-135` | Searched 104 tracked Markdown files and checked 134 unique output JSON/CSV references. | Low |

## J. Other trust risks

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| Two old numbered Python scripts contain unresolved merge-conflict markers and fail syntax parsing. | `src/49_interactive_3d_dashboard.py:390-394`; `src/50_clinical_cdss_app.py:1044-1048` | Python AST parse over 162 in-scope scripts: 160 parse, two fail. | High |
| Legacy files retain unsupported “first”, “clinically viable”, and validation claims despite archive banners. | `archive_old_gbm/most_revolutionary_track.md:9-31`; `archive_old_gbm/OLD_README.md:201,227-247,309-315` | Compared historical claims with current status and later ledger corrections. | High |

## K. Old work claims

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| Pre-2br ledger and old Markdown contain retracted or qualified claims: random-cell C-GAT leakage; assumed kill rates; no demonstrated PDE gain over geometric baseline; stale ρ cutoff. | `vault/PAPER_FINDINGS_LEDGER.md:11,20,31,77`; `vault/README_STALE_CLAIMS.md:17`; `README_STATUS.md:39,46` | Checked representative high-risk entries and cited outputs; complete claim-by-claim check not achieved. See old-tracks note. | High |
| Archived Track C docs still overclaim firstness, clinical viability, and validation. | `archive_old_gbm/most_revolutionary_track.md:9-31`; `archive_old_gbm/OLD_README.md:201,227-247,309-315` | Later status files contradict these framings and describe simulation/assumed parameters. | High |

## L. Old code

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| Two scripts fail syntax parse; 160 of 162 inspected scripts parse. | `src/49_interactive_3d_dashboard.py:390-394`; `src/50_clinical_cdss_app.py:1044-1048` | AST parsing only; runtime imports not tested. | High |
| Older scripts use working-directory-relative output paths against `CLAUDE.md` rule. Examples include scripts 15, 18, 20, 49, and 58. | `CLAUDE.md:81-94`; `src/15_baseline_scvi.py:40-41`; `src/58_rl_adaptive_steering.py:71-72` | Static source inspection; no heavy jobs run. | Medium |
| `output/h3v_test_pairs.csv` and `output/drug_gating_report.md` are named by source but absent. | `src/130_h3v_frozen_lumiere.py:11,29-30`; `src/34_drug_gating_report.py:18-19` | Compared declared paths with filesystem. Script 130’s cited JSON exists. | Low |

## M. Cross-document consistency

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| A21 power loss differs between review overview and saved JSON; the JSON contains `-0.703`, while overview stops at `-0.46`. | `paper/overview_for_review_rano.md:29`; `vault/PAPER_FINDINGS_LEDGER.md:130`; `output/conversion_audit/conversion_audit_summary.json` | Direct value comparison. | High |
| A15 and A18 report upper CI limits 63.8% and 63.5% from separate bootstrap results, not one shared interval. | `README.md:11`; `output/rano_noise/rano_tool_dependence.json`; `output/rano_sensitivity/rano_sensitivity_summary.json` | A15 upper endpoint 0.638095; A18 floor-500 endpoint 0.635419. | Low |

## N. Dead files and large files

| Finding | File:line | Evidence | Severity |
|---|---|---|---|
| Tracked data file exceeds 10 MB: 24.42 MB. | `data/tcga_gbm_expression.tsv` | `git ls-files` and filesystem size. | Low |
| Untracked caches, logs, run packs, viewer output, and `paper/_apply7.py` appear abandoned or generated. | `git status --short`; `paper/_apply7.py` | Worktree inventory; no cleanup performed. | Low |
| 17 old output references in archived plan have no corresponding files. | `archive_old_gbm/paper/GRAND_PLAN.md:90-135` | Reference sweep listed in I and docs/output note. | Low |

## Cannot verify

- Reproducibility of scripts 140–145 to three decimals; no in-place runs were made.
- Runtime importability or output reproduction of old scripts; heavy scripts were not run.
- Every old Markdown claim and every pre-2br ledger citation; K review is high-risk sampling, not exhaustive.
- Whether any person accessed H-2/H-3 outcomes manually before freeze; repository evidence cannot prove that history.
- Clinical validity, generalization, and the true measurement-error distribution from the assumed simulation and single repeatability study.

## Questions for the author

- Should the seven A21 simulation values be described as a selected grid, or should a future declared analysis cover all 18 conversion values?
- Should `START_HERE.md`, the plan, overview, and problem brief be updated to mark calculator and viewer complete and correct the A21 power range?
- Should the two syntax-broken historical scripts remain in the archive, or be excluded from any runnable script inventory?
