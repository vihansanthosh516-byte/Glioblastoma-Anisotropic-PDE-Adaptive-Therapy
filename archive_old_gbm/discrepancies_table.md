> **STATUS 2026-10-18: HISTORICAL.** Written before the pivot to the probabilistic-RANO project. Goals, plans, and words like "first", "novel", "revolutionary", "validated" in this file are NOT current and NOT supported. Numbers are records of past runs; check `vault/PAPER_FINDINGS_LEDGER.md` before you quote one. Current state: `START_HERE.md`.

# Discrepancies Between ISEF Outline and Actual Repo Data

This table checks 12 selected historical outline claims; the exact outline version was not cited, so completeness and original wording are unverified. Corrected 2026-10-19 by reading the full named JSON artifacts. YES/PARTIAL/NO describe historical numerical agreement, not biological or clinical validation. Root output paths below are relative to the repository, not this archive.

| Claim in outline | Source file | Actual value | Match? | Fix |
|---|---|---|---|---|
| **Best single KO SDE2 (C=0.0198)** | output/single_ko_results.json | Full 200-row JSON includes SDE2 collapse_score=0.0; highest simulated single-KO by collapse_score is TMLHE (0.0408) | NO | SDE2 collapse_score is 0.0, not 0.0198; best single KO is TMLHE (collapse_score=0.0408), not SDE2 |
| **Best dual S100A11+ZNF106** | output/dual_ko_results.json | S100A11+ZNF106 pair NOT found; first listed pair: S100A11+S100A6 (C_A=0.003971, C_B=0.003374); second S100A6+MT-ATP6 | NO | S100A11+ZNF106 not present in dual_ko_results.json; First listed pair is S100A11+S100A6; all 21 rows checked, none contains ZNF106. These are model scores, not drug efficacy. |
| **373 directed edges** | output/grn_metrics.json | n_edges: 320 | NO | Actual is 320 directed edges, not 373 |
| **Master switches APOD(46), S100B(45), MT3(40)** | output/grn_metrics.json, output/master_switches.tsv | APOD: n_targets=46 ✓; S100B: n_targets=42 ✗ (outline says 45); MT3: n_targets=40 ✓ | NO | S100B out_degree/n_targets is 42, not 45. APOD(46) and MT3(40) are correct |
| **RL 1.04 mm³ vs Stupp 11.01 mm³** | output/phase5_adaptive_metrics.json | RL final_volume_mm3=0.247, Stupp final_volume_mm3=11.006 | NO | Actual simulated RL final volume is 0.247 mm³ (not 1.04); Stupp 11.01 matches. Assumed treatment model, no measured patient outcome. |
| **C-index for TCGA expression cohort** | output/penalized_survival_metrics.json | c_index: 0.6385432126955132 | PARTIAL | Descriptive internal C-index 0.6385 (n=150); outline expectation/split not supplied. No gene-coefficient p-values in this JSON. |
| **TCGA survival n=518, age HR, p-value** | vault/PAPER_FINDINGS_LEDGER.md A3/2; output/tcga_km_summary.json | n=518 patients, age HR=1.031 (95% CI 1.023-1.038), p=8.88e-16 | YES | Historical age association matches ledger A3; not independent validation. TCGA analyzed cohort: 518, 428 events, 90 censored; KM median 432 days [393,457]. CGGA pooled age HR 1.009, Holm p=0.082 (ledger 2k). |
| **Subtype OS values** | output/subtype_recurrence_summary.json | Proneural: 405.0 days, Classical: 388.0 days, Mesenchymal: 338.5 days, Neural: 261.0 days | PARTIAL | Outline claims to verify but does not specify expected values; actual values reported above |
| **Dose-response TI values** | output/dose_response_report.json | S100A6: TI=5.44, S100A11: TI=5.24, S100A8: TI=4.66, CCL3L1: TI=3.72 | PARTIAL | Outline claims to verify but does not specify expected values; Formula-derived model selectivity scores reported above; not clinical therapeutic indices or efficacy gates |
| **CSGT p-value** | output/csgt_metrics.json | kruskal_wallis p_value: 1.6846997820703238e-31 | YES | Outline says p<0.001; actual 1.68e-31 < 0.001; H=141.7171. Group differences do not prove a continuous transition |
| **SPIB saddle validation** | output/spib_saddle_point_metrics.json | overall_pass: true | YES | Internal numerical consistency flag is true in a synthetic latent landscape; not proof of a biological saddle. |
| **C-GAT accuracy 78.7%** | output/cgat/gat_metrics.json | accuracy: 0.7873333333333333 = 78.73% | YES | Original leaky random-cell score only; patient-level leak-free accuracy 63.8% has donor/composition confounding (ledger 2o). Not a biology result. |

## Historical comparison counts

Twelve claims checked: four YES (33.3%), three PARTIAL (25.0%), five NO (41.7%).
YES plus PARTIAL totals seven (58.3%), not seven full matches or 50%.
These are counts of the retained comparison labels; missing outline provenance
prevents a claim that every factual statement was checked or verified.

Subtype values are descriptive recorded medians; Mesenchymal is 338.5 days,
not 338 under an unspecified rounding rule. Simulated knockout, RL, inferred
network and selectivity outputs do not establish treatment efficacy. The
separate TCGA model is observational, and its internal C-index does not validate
clinical use. See `../START_HERE.md` for the current RANO goal.
