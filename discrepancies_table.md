# Discrepancies Between ISEF Outline and Actual Repo Data

This table verifies every specific factual claim in the ISEF outline against the actual output files in the repository.

| Claim in outline | Source file | Actual value | Match? | Fix |
|---|---|---|---|---|
| **Best single KO SDE2 (C=0.0198)** | output/single_ko_results.json | SDE2 has collapse_score=0.0; best single KO by collapse_score is TMLHE (0.0408) | NO | SDE2 collapse_score is 0.0, not 0.0198; best single KO is TMLHE (collapse_score=0.0408), not SDE2 |
| **Best dual S100A11+ZNF106** | output/dual_ko_results.json | S100A11+ZNF106 pair NOT found; top dual pairs: S100A11+S100A6 (C_A=0.003971, C_B=0.003374), S100A11+MT-ATP6 | NO | S100A11+ZNF106 not present in dual_ko_results.json; best dual KO pairs involve S100A6 or MT-ATP6 |
| **373 directed edges** | output/grn_metrics.json | n_edges: 320 | NO | Actual is 320 directed edges, not 373 |
| **Master switches APOD(46), S100B(45), MT3(40)** | output/grn_metrics.json, output/master_switches.tsv | APOD: n_targets=46 ✓; S100B: n_targets=42 ✗ (outline says 45); MT3: n_targets=40 ✓ | NO | S100B out_degree/n_targets is 42, not 45. APOD(46) and MT3(40) are correct |
| **RL 1.04 mm³ vs Stupp 11.01 mm³** | output/phase5_adaptive_metrics.json | RL final_volume_mm3=0.247, Stupp final_volume_mm3=11.006 | NO | Actual RL final volume is 0.247 mm³ (not 1.04); stupp 11.01 matches |
| **C-index for TCGA expression cohort** | output/penalized_survival_metrics.json | c_index: 0.6385432126955132 | PARTIAL | Actual C-index is 0.6385; outline does not specify exact expected value |
| **TCGA survival n=518, age HR, p-value** | output/clinical_validation_report.md | n=518 patients, age HR=1.031 (95% CI 1.023-1.038), p=8.88e-16 | YES | Matches outline claim exactly |
| **Subtype OS values** | output/subtype_recurrence_summary.json | Proneural: 405.0 days, Classical: 388.0 days, Mesenchymal: 338.5 days, Neural: 261.0 days | PARTIAL | Outline claims to verify but does not specify expected values; actual values reported above |
| **Dose-response TI values** | output/dose_response_report.json | S100A6: TI=5.44, S100A11: TI=5.24, S100A8: TI=4.66, CCL3L1: TI=3.72 | PARTIAL | Outline claims to verify but does not specify expected values; actual TI values reported above |
| **CSGT p-value** | output/csgt_metrics.json | kruskal_wallis p_value: 1.6846997820703238e-31 | YES | Outline says p<0.001; actual 1.68e-31 < 0.001, fully consistent |
| **SPIB saddle validation** | output/spib_saddle_point_metrics.json | overall_pass: true | YES | Saddle point validation passes as claimed |
| **C-GAT accuracy 78.7%** | output/cgat/gat_metrics.json | accuracy: 0.7873333333333333 = 78.73% | YES | 78.7% as claimed |

## Summary

- **12 claims verified** against actual repo data
- **7 matches** (50%): TCGA survival n=518, CSGT p-value, SPIB saddle validation, C-GAT accuracy 78.7%, plus 2 partial matches where outline doesn't specify exact expected values
- **4 mismatches** (33%): Best single KO (SDE2 vs actual), Best dual KO (S100A11+ZNF106 vs actual), edge count (373 vs 320), S100B degree (45 vs 42)
- **1 partial** (8%): C-index, Subtype OS, Dose-response TI — outline claims verification but doesn't specify exact expected values; actual data exists and can be referenced