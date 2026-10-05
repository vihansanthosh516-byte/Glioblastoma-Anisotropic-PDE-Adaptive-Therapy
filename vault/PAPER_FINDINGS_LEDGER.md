# Paper: Findings Ledger

The single source for every number in the paper. Written 2026-10-04.
"Verified" means I read the value from the JSON or CSV on disk in this audit. "Vault only" means the number comes from a vault note and I did not re-read the file.

## Main-text findings (top 3 per track)

### Track A
| ID | Claim | Value | Status | Source file | Must also say |
|---|---|---|---|---|---|
| A1 | C-GAT beats scVI at zone classification | 78.7% vs 73.0% accuracy (+5.8 pp); macro F1 0.783 vs 0.729; AUC 0.922 vs 0.880 | **Verified** (C-GAT and scVI from JSON) | `output/cgat/gat_metrics.json`, `benchmark_comparison.tsv` | 15,000-cell subsample, 3 classes. One seed, no CI. Script 14 was skipped, so the leaderboard was merged by hand. RF and LR are 72.6% and 69.8%. |
| A2 | SPIB index-1 saddle, 5/5 checks | saddle energy 2.019; Hessian eigenvalues +2/-1; attractor energies 0.016 (Healthy), 0.008 (Core), 0.0 (Periphery) | **Verified** | `output/spib_saddle_point_metrics.json` | SPIB found only 2 metastable states (Core+Periphery merged). The saddle is between that merged basin and Healthy. "5/5" are internal consistency checks, not a test on new data. |
| A3 | Age is the dominant prognostic factor in TCGA-GBM | n=518, HR 1.031 (1.023-1.038), p = 8.88e-16 | **Verified** | `output/survival_stats_summary.json`, `tcga_km_summary.json` | Gender and subtype are null (HR 1.15, p 0.17; 1.09, p 0.50). Age is a known factor, so this validates the pipeline; it is not a discovery. Cohort: 428 events, 90 censored; KM median OS 432 d (393-457) (`output/tcga_km_summary.json`). |

Supplement-only Track A: CSGT (H = 141.717, p = 1.68e-31 from JSON; note Periphery T-score 0.649 is not below Core 0.647), GRN (320 edges; master switches APOD 46, S100B 42, MT3 40; 37/380 edges survive bootstrap), three invasion models (26.6, 18.1, 26.6 um/hr), drug screen (MT-CO2 TI +4.82; no dual pair reaches TI > 10), penalized Cox C-index 0.639 (n=150).

### Track B
| ID | Claim | Value | Status | Source file | Must also say |
|---|---|---|---|---|---|
| B1 | Adaptive therapy selects less resistance than MTD | 54/61 patients; 0/61 reverse; median resistant fraction 99.2% (MTD) vs 30.6% (adaptive); adaptive uses 32.4% of MTD dose | **Verified** (n, 32.36%, 10/61, 39/61 from JSON; 54/61 and 99.2/30.6 vault only) | `output/adaptive_cohort_summary.json` | Adaptive ends with **higher final tumour mass** (315.9 vs 306.4; VR 2.39). 10/61 progressed earlier. 39/61 censored in both arms. Resistance in this model is a modelling assumption (E_MAX = rho x 1000). 42 negative-rho patients were excluded. |
| B2 | Forecast beats no-change on the cellular core | Dice 0.259 vs 0.233; diff +0.026 (CI 0.008-0.044); Wilcoxon p 0.0007, **Holm p 0.0015**; n=133 | **Verified** | `output/forecast_labels_core/results.json` | Median diff is only +0.005. Better in 55% of patients. Loses on shrinking cores (-0.039, n=50). DTI adds nothing (-0.001, p=0.97). Target was chosen after a volumes-only probe (disclosed in script header). 19 patients skipped (empty core). The earlier `mask>0` target was negative. Absolute Dice is low. |
| B3 | Anisotropy elongates the simulated tumour | 1.21 vs 1.03; diff +0.178; 98/103; dz 2.89; p = 4.1e-19 (150 d) | **Verified** | `output/fractal_aniso_vs_iso.json` | It is a simulation-vs-simulation comparison, not a test against real tumours. D_f goes the wrong way. Tract alignment is lower under aniso (0.45 vs 0.50). |

Supplement-only Track B: Sobol (PDE: alpha_sens 0.43, rho 0.20, D_w 0.000; N=128), MPC (8 patients), 3D (8 patients, sparing 81.4% +- 30.4%), horizon crossover (Dice 0.923 at D=0.1, 90 d), script 77/80 ablation.

### Track C
| ID | Claim | Value | Status | Source file | Must also say |
|---|---|---|---|---|---|
| C1 | Kill-rate-conditioned policy wins at equal budget | 100% win (cohort64, lhs60, real_test); blind PPO 25% on cohort64 | **Verified** | `output/rl_kill_conditioned/evaluate.json` | Win is on **day-90 volume**, which never wins on TTP (0.0 / -29.4 / -10.7 / -7.9 d across the four sets, script 68). A one-line efficiency rule matches the RL policy (0.2% of oracle). So RL adds nothing over the rule. |
| C2 | Probe-then-commit | 100% noise-free; 62% (sigma 0.10, rule arm, 3-seed mean 45/75/65) | **Verified** | `output/probe_paced_policies/evaluate.json` | cohort64 only. Falls to 11-46% on other sets under noise (vault only). Identifiability result, not a policy. |
| C3 | Pacing extends time to progression | TTP +48.5 / +52.5 / +2.1 / +15.7 d; longer in 101/111, shorter in 0 | **Verified** | `output/probe_paced_policies/evaluate.json` | **Real-parameter set gains only +2.1 d.** Day-90 volume loses. It is pacing, not intelligence. Report all four sets. |

Supplement-only Track C: inverse estimation (124/154 at lower bound), robust MPC, ablation, biomarker, virtual cohort, plus the three resistance negatives (scripts 76, 82, 83).

## Limitations section (from the vault, verified)
| Limitation | Number | Source |
|---|---|---|
| Shrinkage cannot be forecast | 80/152 shrank or stayed (53%) | `forecast_validation/forecast_results.json` via `NEGATIVES_REVISITED.md` |
| Forecast loses on shrinking cores | -0.039 Dice (n=50) | `forecast_labels_core/results.json` |
| MGMT does not predict TTP | 169 vs 181 d; log-rank p = 0.5998; n=130 | `output/mgmt_resistance/results.json` (**Verified**) |
| Resistance-driven adaptive gain fails on real patients | 0 of 6 cells (script 83); no CI above 0 on real_test (script 82); 2 of 384 rows (script 76) | `plastic_resistance/results.json`, `mgmt_resistance/results.json` |
| D_f is not valid | wrong direction | `fractal_aniso_vs_iso.json` |
| Script 60 used D about 10x too low | 0.013 vs 0.13 mm2/d | [[Script-78-Horizon-Crossover]] |
| 124/154 patients fit at the growth lower bound | | [[TRACK_C_RESULTS]] |

## Open conflicts (fix before the paper quotes them)
1. ~~TCGA events and median OS~~ **RESOLVED 2026-10-04 (script 84).** n=518, 428 events, 90 censored, **KM median 432 d (95% CI 393-457)**. Source: `output/tcga_km_summary.json`. The old 377.5 d (all times) and 383 d (deaths only) are not survival medians. The raw TCGA CSV flags 441 deaths; the analysed cohort file (`clinical_mapped_cohort.csv`) has 428. State which file the paper uses (the cohort file).
2. ~~Biomarker threshold~~ **RESOLVED.** Script 62 is current: rho* = 0.0751 /day, bootstrap CI [0.0746, 0.0820], early start wins for 60/64, only 4 patients above rho*. README 0.024 is an older analysis.
2b. ~~Decisions.md "23 of 61"~~ **RESOLVED (script 85).** 10 patients have rho > 0.02; all 10 progress earlier under adaptive; 0 of 51 below. Usable as a stratified B1 finding, with the caveats that n=10 and the cutoff was chosen after seeing the data. Source: `output/rho_threshold_check.json`.
2c. **S100A6 therapeutic index (5.44, highest) uses a Cox weight of about -2e-7 (zero).** Do not claim S100A6 is prognostic. Only S100A8 (0.043) and CCL3L1 (0.021) have non-zero directional weights; none is significant.
2d. **Script 47 dual-agent uses 2.5-3.3x more drug than MTD** (AUC 944-1280 vs 385). Do not call it dose-sparing. Single-agent is earlier than MTD for 1/8 patients, not 4/8.
2f. **Re-checked 2026-10-04 and found correct:** Script 27 necrosis snapshots (22.2 / 30.4 / 38.1%); dual-KO screen (21 pairs, no positive Bliss, TI range -0.98 to +0.44); Script 76 (19 of 384 rows beat best non-adaptive; 17 are cohort64; 2 lhs60 rows +22.3 and +15.9 d); Script 78 Dice table (all 16 values for rho 0.02 at 2 mm and the 1 mm checks).
2g. **Not reproducible, do not cite:** Script 27 "front velocity 0.118-0.174".
2h. **Still not re-checked:** Script 29 hand-computed free-propagation numbers beyond `invasion_summary.json` (26.6 um/hr and 28.3 analytical are verified; the 5.97% error was not recomputed).
2e. **Script 48 (3D) uses legacy D** (0.013 / 0.0013 mm2/d), about 10x below literature.
3. **"97 scripts".** True count of numbered files in `src/` is 97, but they use 82 distinct numbers (max 83). Say "97 numbered scripts".
4. **Script 66 numbers** (PPO 12.67 vs Stupp 14.47 mm3) differ from README (13.94 vs 11.01). Use script 66.
5. **Track A numbers I did not re-read:** scVI 0.731, NMF, GRN bootstrap (37/380), drug-screen TI values, penalized Cox. Treat as vault-only.
