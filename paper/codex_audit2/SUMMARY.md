# Round 2 audit summary

All 74 fixed chunks are DONE. Chunk 004 was split into 004a and 004b; both reports exist. The original chunk_004.md is retained as a superseded draft and is excluded from the claim counts below to avoid double-counting. The 76 progress rows include the two split reports. No PARTIAL chunks.

The audit was static except for `py_compile` and the required pytest run in a scratch copy. No project script was executed. Status counts below count claim-table rows in the final reports; script rows use the requested severity scale instead of claim labels. A claim can have more than one issue row.

## Counts

| Kind / track bucket | CORRECT | WRONG_NUMBER | NO_SOURCE | STALE | OVERCLAIM | LIMIT_MISSING | UNVERIFIABLE | High | Medium | Low |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Ledger, old entries 1–2bq (chunks 001–013, 004a/b) | 400 | 18 | 16 | 39 | 13 | 39 | 68 | — | — | — |
| Markdown, archived old-GBM files (chunks 017–033; includes mixed review-response file) | 104 | 56 | 17 | 97 | 74 | 123 | 153 | — | — | — |
| Markdown, current/background/audit files (chunks 014–016, 034–045) | 195 | 19 | 2 | 89 | 7 | 111 | 151 | — | — | — |
| Markdown total (chunks 014–045) | 299 | 75 | 19 | 186 | 81 | 234 | 304 | — | — | — |
| Scripts, Track A (numbered 01–41) | — | — | — | — | — | — | — | 56 | 54 | 68 |
| Scripts, Track B (42–48) | — | — | — | — | — | — | — | 7 | 9 | 5 |
| Scripts, Track C/C-Extended (49–65) | — | — | — | — | — | — | — | 24 | 29 | 32 |
| Scripts, forecast/PredictGBM (66–139) | — | — | — | — | — | — | — | 19 | 82 | 152 |
| Scripts, current RANO (140–145) | — | — | — | — | — | — | — | 0 | 11 | 9 |
| Scripts, utility/unclassified paths | — | — | — | — | — | — | — | 37 | 49 | 51 |
| Tests (chunk 074) | Read 17 test files; scratch `pytest tests/ -q`: 108 passed, 13 warnings, exit 0 | | | | | | | | | |

Track assignment for scripts uses the numbered script prefix in each report. “Utility/unclassified” covers unnumbered files and rows without an unambiguous track number. Archived Markdown track buckets are based on the fixed chunk file set; several files discuss more than one track. The current/background bucket is not a scientific track.

## 30 worst findings

1. **Cell-level leakage in Track A classifiers.** `src/04_export_for_attention_model.py:102–112` makes a stratified random cell split; patient identity crosses train/test. The ledger notes patient lookup reaches 82.8%. Evidence: chunk_046, `output/patient_level_cv.json`, ledger A1.
2. **Repeated test-set model selection.** Scripts 05, 07, 08 inspect test metrics every epoch and report best test accuracy. Evidence: chunk_046, `src/05_attention_gated_network.py:153–184`, `src/07_method2_transformer.py:142–174`, `src/08_method3_hybrid.py:205–239`.
3. **GAT test data used for checkpoint choice.** `src/13_gat_train.py:150–185` selects the checkpoint by test accuracy. Evidence: chunk_052.
4. **Graph includes label-derived region features.** `src/12_gat_build_graph.py:36–47,65–97` uses tissue histology to encode `same_region`, transition type and region difference. Evidence: chunk_051.
5. **Reported model size wrong.** Script 05 says ~3.2M parameters; `output/nn_metrics.json:n_params` is 488,899. Evidence: chunk_046.
6. **Transformer feature count wrong.** Script 07 says 500 HVGs; code sets `N_GENES=100`, and selected columns are not ranked by HVG importance. Evidence: chunk_046.
7. **Benchmark script should not run as a current analysis.** `src/09_benchmark_comparison.py` writes relative to cwd, recreates an archived report, emits unsupported “novel” language and reports test-selected best-epoch metrics. Evidence: chunk_047.
8. **Wrong Track A comparison label.** Script 02’s `Core_vs_Peri` output is Core-vs-rest; script 03 propagates the wrong label and hard-codes marker lists. Evidence: chunk_046, `src/02_preprocess_umap_de.py:207–224`.
9. **A6.2 bootstrap count mismatch.** Script 103 uses 2,000 bootstrap draws although the later plan amendment requires 10,000. Evidence: chunk_047, output commit date and `analysis_plan_v1.md` A6.2.
10. **H-2 test freeze is not one-shot enforced.** Script 123 scores the spent set without one-shot enforcement. Evidence: chunk_050; H-2 is marked spent in `analysis_plan_v1.md` and `START_HERE.md`.
11. **H-3 confirmation set overlaps model selection.** Script 125’s “8-patient confirmation” includes two patients used to select the variant. Evidence: chunk_050.
12. **H-3v sample count wrong and config may be overwritten.** Script 127 says n=105; output/ledger show n=104. Its freeze writer can overwrite config. Evidence: chunk_050.
13. **Final H-3 arm ranking contradicts the stated rule.** Script 132 selects `hybrid_mri_rank_mu` while listed `hybrid_rank_mu` has the higher OOF value (0.8529 vs 0.8519). Evidence: chunk_051, `final_train.json`.
14. **H-3 locked score lacks full input/run protections.** Script 133 does not enforce `runs_allowed: 1` or hash ALPTDG inputs/context logits. Evidence: chunk_051.
15. **Script 95 identifies graph/cVAE label leakage.** Its own audit notes all-cell graph/cVAE features can expose held-out labels. Evidence: chunk_065.
16. **Script 31 loads a `.npy` file using `torch.load`.** This is a likely runtime failure / incorrect loader path. Evidence: chunk_056, `src/31_virtual_knockout_engine.py`.
17. **Script 35 synthesizes expression values against TCGA clinical data.** The resulting association is not an observed matched-expression analysis. Evidence: chunk_056.
18. **Script 37 overstates clinical validation.** Its output interpretation is hard-coded as “clinical validation” without a validating cohort or clinical endpoint. Evidence: chunk_056.
19. **Script 34 carries stale typed result constants.** Printed/reported values are not computed from the current saved outputs. Evidence: chunk_056.
20. **Several Track C results are synthetic but described clinically.** Script 59 and neighboring solvers use hard-coded PDE/therapy constants; expected outputs for scripts 55–57 are absent. Evidence: chunks 059–060.
21. **HIL latency claim is not guaranteed.** `src/hil/pump_interface.py` claims “<500ms guarantee”; code does not enforce a timing bound, and serial dependency is undeclared. Evidence: chunk_068.
22. **FNO benchmark uses random input and hard-coded baseline timing.** `src/neural_pde/benchmark_fno.py` uses random inputs and types a 3-second FD estimate. Evidence: chunk_068.
23. **FNO training targets are not a PDE solver result.** `src/neural_pde/generate_fno_data.py` uses a simple exponential update. Evidence: chunk_069.
24. **`predictgbm_io.py` references `os` without importing it.** Evidence: chunk_069.
25. **Pseudotrajectories are extrapolated, not observed patient trajectories.** `src/pseudotrajectory_generator.py` should not be described as longitudinal clinical data. Evidence: chunk_069.
26. **Exploratory post-freeze analysis reads spent H-2 test outputs.** `src/run_sensitivity_shifts.py` is explicitly post-freeze exploratory and writes beside those outputs. Exclude from locked-run inventory. Evidence: chunk_071.
27. **SPIB split does not preserve an untouched test set.** Training receives the test dataset in `src/spib_train.py`. Evidence: chunk_072.
28. **Treatment infusion simulation is described causally without clinical causal evidence.** `src/treatment_aware_pde.py` is a model simulation. Evidence: chunk_072.
29. **Script 09/empty batch runner/graph builder have unsafe or absent implementation.** `src/batch_run_digital_twins.py` is 0 bytes; `src/build_full_graph.py` does heavy work at import; script 09 writes legacy outputs. Evidence: chunks 047, 066.
30. **The 2026-10-04 status table in `vault/LIMITS_FIX_OPTIONS.md` is stale.** It says leak-free C-GAT is “Running,” while `output/cgat_leak_free.json` and ledger A1 hold completed results; L24’s conditional is superseded by all 61 scores being 1.0/no A-to-B link. Evidence: chunk_036, ledger A24.

## Ledger entries needing correction

The full row-by-row evidence and one-line correction is in each cited chunk report. Every non-CORRECT entry is listed here; entries with multiple claim rows have all claim-specific fixes in that report.

| Entry | Correction / report |
|---|---|
| 1 | Cite raw CSV path/event field for 441 deaths, or remove it. (chunk_001) |
| 2 | Identify README source for 0.024 and label historical, or remove. (001) |
| 2b | Replace 23/61 with script-85 stratification and mark superseded. (001) |
| 2c | Remove S100A6 5.44/highest unless dated source is identified; current maximum is MT-CO2 4.8193. (001) |
| 2e | Cite comparable literature parameters or remove “below literature.” (001) |
| 2f | Correct count to 20 snapshot values; identify Dice comparison. (001) |
| 2g | Identify source for velocity range or delete it. (002) |
| 2h | Correct relative error wording/value to 5.98%. (002) |
| 2i | Call 358.9 days a horizon-censored simulation mean, not uncensored TTP. (002) |
| 2j | Label drug inputs fitted-growth/assumed-kill, not observed efficacy. (002) |
| 2k | Cite a run record for date; keep pooled c-index out unless validated. (003) |
| 2m | State likelihood qualifier and unconfirmed provenance. (003) |
| 2o | Say the classifier beats this coarse permutation null; do not infer biological signal. (004a) |
| 2q | Label V7b as post hoc sensitivity. (004b) |
| 2r | Replace interim intervals with final 10,000-draw values. (004b) |
| 2s | Say results are compatible with poor persistence under assumed boundary perturbation. (004b) |
| 2u | Replace 93 with 101/101 and keep split-run caveat. (005) |
| 2v | Mark superseded by 2aa; retain only as interim history. (005) |
| 2w | Mark explicitly superseded by 2y. (005) |
| 2x | Add manifest and pair-table row references. (005) |
| 2y | Cite run logs/manifest for exact modes. (005) |
| 2z | Cite declaration/config for stratification and grid. (005) |
| 2aa | State what 120 counts; current output has seven top-level strata. (006) |
| 2ab | Describe as two-tool disagreement, not measured scan-rescan noise. (006) |
| 2ac | Mark NNLS config superseded by 2ae. (006) |
| 2ad | Label result under assumed perturbation, not measured noise. (006) |
| 2ae | Keep exploratory; add scratch log/artifact or omit performance claim. (006) |
| 2af | Scope bound to exploratory volume-matched GBM setup and CI. (006) |
| 2ag | Correct deltas to 1.03→2.46; 0.878→0.876; 0.353→0.353. (007) |
| 2ah | Cite threshold and horizon method. (007) |
| 2ai | Cite dated access record or label author-reported. (007) |
| 2aj | Add exact README path/quote or remove. (007) |
| 2ak | Cite manifest counts. (007) |
| 2al | Cite cache/manifest exclusions. (007) |
| 2am | Retain 0.553; source the range or say unverified. (008) |
| 2an | Cite patient-level skip audit or omit. (008) |
| 2ao | Cite training log; avoid unsupported underfit diagnosis. (008) |
| 2ap | Add input inventory or label synthetic fallback. (008) |
| 2aq | Cite code implementing barrier threshold. (008) |
| 2ar | Cite meta-analysis; clarify observational OS, not PFS treatment effect. (008) |
| 2as | Cite runtime record or remove. (009) |
| 2at | State tested grid did not fix shrink subgroup; required rates unknown. (009) |
| 2au | Say M1 lowered this development ensemble’s equal-volume mean. (009) |
| 2av | Cite runtime record or remove. (009) |
| 2aw | Keep 15 epochs; state ~446–447 training examples per fold, 476 only for pool. (009) |
| 2ax | Cite freeze/tag and run log for once-only claims. (009) |
| 2az | Change “all” to “178 of 180 have parseable surgery and RT dates.” (010) |
| 2bb | Say “not observed with this proxy in TUM dev.” (010) |
| 2bc | Cite date-anonymization documentation for preserved intervals. (010) |
| 2bd | Cite manifest keys for seed, sample, epochs, runtime, hashes. (010) |
| 2be | Say “post hoc pooled development pick; requires new locked test.” (011) |
| 2bf | Keep as dated preparation record; point to 2bh for outcome. (011) |
| 2bg | Attach attempt-1 Kaggle log/export or mark author-reported. (011) |
| 2bi | Attach comparison output or say artifact unavailable. (011) |
| 2bj | Recount from explicit scope; update status docs. (011) |
| 2bk | Cite reading log or call review author-reported. (012) |
| 2bl | Replace expectation with measured threshold sweep; mark superseded by 2bn. (012) |
| 2bm | Keep as author-reported process or cite dated access log. (012) |
| 2bn | Correct SHA-256 to `66d7b40b9cb76c6203a8ea8b1cad89cc67537951533796a7fa2d2568cde713d8`. (012) |
| 2bo | Label as self-reported run provenance or retain scratch manifest/transcript. (012) |
| 2bp | Correct HTML size to 15.8 MB or state units/date. (012) |
| 2bq | Keep “run once” only with execution log. (013) |

Entries 2ay, 2ba, 2bh, 2bu–2bx have CORRECT claims in their audited rows; see reports for mixed entries with some correct and some issues. “No correction” here does not mean all claims in an entry are independently validated.

## Markdown files needing correction

Every file below has at least one non-CORRECT claim or a file-level correction. The exact evidence and one-line edits are in the named chunk report; corrections are claim-specific there. Archived files should retain their HISTORICAL banner and unsupported old claims must not be presented as current.

| File | Correction locator |
|---|---|
| `AGENTS.md` | Apply procedural-reference corrections in chunk_014. |
| `CLAUDE.md` | Update stale script/Amendment counts and current command references (014). |
| `README.md` | Correct current summary claims listed in 014. |
| `README_STATUS.md` | Correct blanket status claims and align to `START_HERE.md` (014). |
| `START_HERE.md` | Correct stale inventory/status statements in chunk_015. |
| `analysis_plan_v1.md` | Correct declaration/amendment statements identified in chunk_016. |
| `archive_old_gbm/FIXES_SUMMARY.md` | Use chunk_017 claim-level fixes; retain historical banner. |
| `archive_old_gbm/OLD_README.md` | Use chunk_018 claim-level fixes; retain historical banner. |
| `archive_old_gbm/README.md` | Use chunk_018 claim-level fixes; retain historical banner. |
| `archive_old_gbm/REPO_MIND_MAP.md` | Correct unsupported relationships/claims in chunk_019; retain banner. |
| `archive_old_gbm/TRACK_A_RESULTS.md` | Correct claims in chunk_020; retain banner. |
| `archive_old_gbm/discrepancies_table.md` | Correct claims in chunk_020; retain banner. |
| `archive_old_gbm/docs/PROJECT_REVIEW_SUMMARY.md` | Correct claims in chunk_021; retain banner. |
| `archive_old_gbm/docs/SCRIPTS_01_18_ANALYSIS.md` | Correct claims in chunk_022; retain banner. |
| `archive_old_gbm/docs/VALIDATION_SUMMARY.md` | Correct count conflicts and missing artifact citations in chunk_023; retain banner. |
| `archive_old_gbm/docs/article_summary.md` | Verify/correct literature details in chunk_023; retain banner. |
| `archive_old_gbm/docs/dataset_info.md` | Correct cohort counts/definitions in chunk_023; retain banner. |
| `archive_old_gbm/docs/isef_compliance.md` | Update stale deadlines/rules and verify cited policy in chunk_023; retain banner. |
| `archive_old_gbm/docs/methodology_upgrade_summary.md` | Correct claims in chunk_024; retain banner. |
| `archive_old_gbm/most_revolutionary_track.md` | Remove unsupported “revolutionary” framing and correct claims (024). |
| `archive_old_gbm/output/BENCHMARK_SUMMARY.md` | Correct benchmark values/limitations in chunk_024; retain banner. |
| `archive_old_gbm/output/EXECUTIVE_SUMMARY.md` | Correct value/provenance claims in chunk_024; retain banner. |
| `archive_old_gbm/output/MONTH10_AUDIT.md` | Correct dated/status claims in chunk_024; retain banner. |
| `archive_old_gbm/output/POSTER_KEY_FINDINGS.md` | Correct claims and remove unsupported clinical framing (025). |
| `archive_old_gbm/output/RESEARCH_EXECUTIVE_SUMMARY_AND_TECHNICAL_DOSSIER.md` | Correct evidence/provenance and validation language (025). |
| `archive_old_gbm/output/clinical_actionability_report.md` | Remove unsupported actionability claims and correct numbers (025). |
| `archive_old_gbm/output/clinical_validation_report.md` | Remove unsupported clinical validation and qualify simulation (025). |
| `archive_old_gbm/output/drug_gating_report.md` | Correct assumptions and clinical wording (025). |
| `archive_old_gbm/paper/GRAND_PLAN.md` | Mark as historical plan; correct stale status/results (026). |
| `archive_old_gbm/paper/RESEARCH_PROMPT_digital_twin.md` | Correct historical claims in chunk_027; retain banner. |
| `archive_old_gbm/paper/masterplan_implementation_plan.md` | Correct plan/status claims in chunk_027; retain banner. |
| `archive_old_gbm/paper/overview_for_review.md` | Correct evidence and results in chunk_027; retain banner. |
| `archive_old_gbm/paper/project_plan_trusted_twin.md` | Correct plan/status claims in chunk_027; retain banner. |
| `archive_old_gbm/paper/review_packet_phase0-3.md` | Correct dated/review claims in chunk_028; retain banner. |
| `archive_old_gbm/paper/summary_all_tracks.md` | Correct stale summary figures and source references (028). |
| `archive_old_gbm/vault/Repo.md` | Correct archive inventory/claims (029); retain banner. |
| `archive_old_gbm/vault/Script-44-Adaptive-Therapy.md` | Correct simulation/assumption claims (029); retain banner. |
| `archive_old_gbm/vault/Script-46-Sensitivity.md` | Correct claims (029); retain banner. |
| `archive_old_gbm/vault/Script-47-Optimal-Control.md` | Correct dose-comparator and efficacy claims (029); retain banner. |
| `archive_old_gbm/vault/Script-48-3D-Extension.md` | Correct diffusivity/validation claims (029); retain banner. |
| `archive_old_gbm/vault/Script-60-66-Swanson-D.md` | Correct provenance and run claims (030); retain banner. |
| `archive_old_gbm/vault/Script-66-RL-Equal-Budget.md` | Correct simulated-vs-real framing (030); retain banner. |
| `archive_old_gbm/vault/Script-67-68-Conditioned-TTP.md` | Correct comparator and TTP claims (030); retain banner. |
| `archive_old_gbm/vault/Script-74-Fractal-Elongation.md` | Correct alignment interpretation (030); retain banner. |
| `archive_old_gbm/vault/Script-75-Probe-Paced.md` | Correct parameter/validation framing (030); retain banner. |
| `archive_old_gbm/vault/Script-78-Horizon-Crossover.md` | Limit conclusion to tested deterministic grid (031); retain banner. |
| `archive_old_gbm/vault/Script-80-Sobol-PDE.md` | Correct scope and parameter assumptions (031); retain banner. |
| `archive_old_gbm/vault/Script-81-Label-Correct-Forecast.md` | Verify shrinking-tumour value and artifact (031); retain banner. |
| `archive_old_gbm/vault/Script-82-MGMT-Resistance.md` | Qualify exploratory spent-set use and smoke-test claim (031); retain banner. |
| `archive_old_gbm/vault/Script-83-Plastic-Resistance.md` | Correct simulation/assumption framing (031); retain banner. |
| `archive_old_gbm/vault/TRACK_A_RESULTS.md` | Correct results in chunk_032; retain banner. |
| `archive_old_gbm/vault/TRACK_B_RESULTS.md` | Correct results in chunk_032; retain banner. |
| `archive_old_gbm/vault/TRACK_C_RESULTS.md` | Correct “Verified Results” overstatement and missing artifacts (033); retain banner. |
| `archive_old_gbm/vault/UCSF-PDGM-Tensor-Study.md` | Verify missing CSVs and qualify conclusions (033); retain banner. |
| `paper/chatgpt_review_response.md` | Update stale dated states and independently verify relayed figures (033). |
| `paper/literature_search_digital_twin.md` | Correct literature claims/search scope in chunk_034. |
| `paper/new_project_options.md` | Update old-goal ranking and cite claims in chunk_034. |
| `paper/overview_for_review_rano.md` | Correct results/limitations per chunk_034. |
| `paper/phd_level_plan_probabilistic_rano.md` | Correct current plan claims per chunk_034. |
| `paper/problem_first_rano_noise.md` | Correct evidence and framing per chunk_034. |
| `paper/research_prompt_v2.md` | Correct goal/status/citation claims per chunk_035. |
| `paper/research_v2_quick_pass.md` | Correct old status and claims per chunk_035. |
| `vault/AUDIT_REPORT.md` | Correct audit scope/findings per chunk_035. |
| `vault/Decisions.md` | Correct supersession/status references per chunk_035. |
| `vault/LIMITS_CURRENT.md` | Correct current-limit statements per chunk_035. |
| `vault/LIMITS_FIX_OPTIONS.md` | Replace stale 2026-10-04 status rows, especially C-GAT and L24 (036). |
| `vault/NEGATIVES_REVISITED.md` | Correct comparison/statistic rows and evidence (037). |
| `vault/Open-Questions.md` | Update resolved questions and cite evidence (037). |
| `vault/PAPER_DATA_SOURCES.md` | Verify/correct source inventories and cohort facts (037). |
| `vault/PAPER_FIGURES.md` | Reconcile 14-vs-13 figures and update stale plan block (037). |
| `vault/PAPER_FINDINGS_LEDGER.md` | Correct the non-CORRECT rows for 2br–2bx in chunk_038; older entries are detailed in chunks 001–013. |
| `vault/PAPER_REFERENCES.md` | Verify bibliography/notes count claims (039). |
| `vault/Project-State.md` | Update stale status and missing evidence links (039). |
| `vault/README_STALE_CLAIMS.md` | Correct stale claim index/source pointers (039). |
| `vault/WEAK_POINTS.md` | Update weakness/status values from current ledger (039). |
| `vault/Welcome.md` | Update links/status guidance (039). |
| `vault/audit_2026-10-18/FIX_PLAN_v1_initial.md` | Mark replaced plan clearly and correct size/status claims (040). |
| `vault/audit_2026-10-18/FIX_PLAN_v2_1_approved.md` | Correct approval/provenance statements where unsupported (040). |
| `vault/audit_2026-10-18/H4_dataset_desk_review.md` | Verify dataset claims and cite primary fields (040). |
| `vault/audit_2026-10-18/INDEX.md` | Correct index counts and artifact references (040). |
| `vault/audit_2026-10-18/P4_threshold_sweep_declaration.md` | Correct declaration/provenance statement if not backed by run record (040). |
| `vault/audit_2026-10-18/P4_threshold_sweep_report.md` | Correct named-patient and run-log claims (041). |
| `vault/audit_2026-10-18/README_STATUS_change_log.md` | Add sources for process/size claims (041). |
| `vault/audit_2026-10-18/kill_rate_grid_declaration.md` | Correct run-once/declaration provenance if no log (041). |
| `vault/audit_2026-10-18/kill_rate_grid_report.md` | Cite literature details and qualify patient examples (041). |
| `vault/audit_2026-10-18/literature_check.md` | Add primary tables for ablation/RHUH claims; keep “not found” scoped (041). |
| `vault/audit_2026-10-18/originals/README_STATUS_original_2026-10-17.md` | Add historical banner; fix stale/spent-lock and numeric claims (042). |
| `vault/audit_2026-10-18/originals/REPO_MIND_MAP_original_2026-10-14.md` | Add historical banner; correct test-count and status claims (043). |
| `vault/audit_2026-10-18/plan_review_record.md` | Cite review/provenance for assertions (044). |
| `vault/audit_2026-10-18/readme_status_audit.md` | Correct item-by-item status against current sources (044). |
| `vault/claim_ledger.md` | Correct claim status and evidence links per chunk_044. |
| `vault/failed_hypotheses.md` | Correct hypothesis outcomes and add source pointers (045). |
| `vault/learned_vs_assumed.md` | Correct assumed/learned labels and cite exact outputs (045). |

## Scripts to exclude from the runnable inventory

- `src/09_benchmark_comparison.py`: archived Track A report generator; cwd-relative writes; unsupported “novel” language; test-selected best-epoch metrics (chunk_047).
- `src/batch_run_digital_twins.py`: zero-byte file; no implementation (chunk_066).
- `src/build_full_graph.py`: performs heavy build work at import and uses cwd assumptions; not safe as a general runnable entry (chunk_066).
- `src/22_saddle_point_proof_backup.py`: alternate legacy solver writes the same cwd-relative output as the primary and can overwrite it (chunk_054).
- `src/run_sensitivity_shifts.py`: explicitly post-freeze exploratory access to spent H-2 outputs; keep out of locked/current run inventory (chunk_071).

Additional old Track A/B/C and forecast scripts are historical background per `START_HERE.md:5–6`; do not describe their existence as current pipeline validation. Several need optional ML/GPU dependencies and many write relative to cwd; see chunks 046–073 before selecting any as runnable.

## Could not verify

- Raw TCGA source count of 441 deaths in ledger entry 1: no raw CSV/key named; mapped output has 428.
- Old README value 0.024 (entry 2); cited README/path absent.
- Many author-reported dates, run-once claims, elapsed times, Kaggle imports, and screenshots lack retained run logs/manifests; see entries 2k, 2as/2av, 2ax, 2bq and chunks 040–041.
- Several cited archived Track C/UCSF CSV/JSON files are absent from current output trees (chunks 032–033).
- Literature claims based on snippet/search-summary review only; primary supplement tables needed for some RHUH, ablation, DTI and dataset details (chunks 023, 034, 039–041).
- Whether all `py_compile` checks reflect runtime importability; compilation does not load dependencies. Missing or undeclared dependencies are listed by script in chunks 046–073.
- Whether lock/freeze hashes and one-shot rules were respected for every historical run where logs are absent (chunks 050–052, 064, 071).
- Patient-level leakage-free performance for old Track A classifiers beyond the saved analyses; C-GAT model weights are missing (chunks 046, 051–053).
- Execution count for scripts labeled “run once” where no run manifest exists (chunks 012, 013, 040–041).

## Questions for the author

1. Should the 0-byte `src/batch_run_digital_twins.py` be deleted or kept as a placeholder?
2. Can you supply the source CSV/key for the raw 441 TCGA deaths and the README behind ledger value 0.024?
3. Are all five scripts in the exclusion list still intended to appear in any top-level runnable inventory?
4. Can you provide the missing archived Track C/UCSF output artifacts and run manifests, or should those claims remain explicitly unverifiable?
5. Should the stale 2026-10-04 `vault/LIMITS_FIX_OPTIONS.md` status table be retained as history with a banner, or removed from current navigation?
6. Which optional ML/GPU dependency manifest is authoritative for runnable inventory: `Requirements.txt`, `pyproject.toml`, or the separate Track-A/WSL setup in `CLAUDE.md`?
