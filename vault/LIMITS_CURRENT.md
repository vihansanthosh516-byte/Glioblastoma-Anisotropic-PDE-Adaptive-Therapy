> **HISTORICAL BACKGROUND (2026-10-19 correction):** October 4 Track A–C limitation snapshot. All treatment/RL summaries concern simulations with assumed kill rates, even where growth inputs are measured. No prospective or clinical validation is established. Current probabilistic-RANO limitations/status: START_HERE.md and ledger 2br–2bx. Unkeyed residual literature/access/status claims below are unverified: exact dated inventories or primary source records are missing.

# Limits: historical Track A–C status (2026-10-04, evening)

Original list for the old GBM Limitations section; not the current RANO limitations. The ideas for fixes are in [[LIMITS_FIX_OPTIONS]]. Numbers are in [[PAPER_FINDINGS_LEDGER]].
Status words: **Fixed** = rerun and the limit is gone. **Checked** = tested, the limit is now measured but still there. **Open** = nothing done. **Finding** = it is a result, not a bug.

## Biggest (change what the paper can claim)
| ID | Limit | Status | What we know now |
|---|---|---|---|
| L0 | Track A uses a random cell-level split, with label-derived graph edges and test-set epoch choice | **Checked; script 88 done, script 95 audit done** | Patient-ID lookup scores 82.8%, above C-GAT 78.7%. Patient-level LR 60.3%, RF 62.0%, scVI-LR 66.2% (script 86). Leak-free C-GAT: fold 0 only (clean PCA graph 77.1%, cVAE without label edges 75.7%, published edges 86.4%). All 5 folds: gat_pca_clean 63.8% (SD 13.3), cVAE no label edges 74.2%, published edges 75.8%. A classifier that sees only cells-per-patient scores 64.6% (SD 27.5). See ledger 2o. |
| L24 | The Track A to Track B link | **Checked: none** | Historical constant-score/ID-mismatch report is unverified here: the exact linked score artifact/key is missing. Track A stands alone. |
| L14/L28 | `real_test` is not a real-drug test | **Checked** | Real growth rates, but the same assumed kill rates for all 21 patients. Median kill/rho 77.7 vs 2.7-12.5 in synthetic sets (script 87). This is a model-regime difference; it does not establish a causal explanation of the +2.1 d gain. |
| L10 | Script 47 dual-agent "rescue" | **Fixed (explained)** | MTD plus the same second drug: 358.9 d, AUC 1152. Dual-adaptive: 305.6 d, AUC 1015 (script 89). The matched second-drug MTD baseline outperforms adaptive in this simulation; it does not establish a patient mechanism. |

## Data
| ID | Limit | Status | Note |
|---|---|---|---|
| L1 | Cellular-core mask forecast: mean Dice 0.2591 vs persistence 0.2332 | Open | Historical selected-target OOF result, `output/forecast_labels_core/results.json` (labels 1+3), not independent validation. Later manifest PDE gain has CI spanning zero (ledger 2w/2y). Old download/size status unverified without dated inventory. |
| L2 | Forecast loses on shrinking tumours (53% shrank or stayed) | **Finding, replicated** | LUMIERE volume forecasting (not mask Dice), `output/lumiere_volume_forecast.json`, ledger 2l: growth-only model loses to no-change (n=40 p 0.0003; n=50 p 0.004; script 91). |
| L3 | UCSF-PDGM: grade 2-3, one scan | Open | RHUH-GBM downloaded (n=40, 3 time points, segmentations), historical status superseded by H-2 scoring and H-4 training (ledger 2ax/2bf); H-2 now spent. |
| L4 | Only 21 real patients in Track C | Open | LUMIERE gives 58-64 growth rates (script 91) but they are not run through Track C. |
| L5 | D not identifiable; 124/154 rho fits at the lower bound (`output/inverse_est_metrics.json:summary`) | Open | |
| L6 | 15,000-cell subsample, one seed, no CI | Partly | Scripts 86/88 add 5-fold mean and SD. |
| L7 | TCGA n=150 with expression; age is known | **Checked** | CGGA GBM (n=374): age HR 1.009 (1.001-1.018), Holm p 0.08, vs TCGA 1.031. Sex, MGMT null (script 90). Gene models not replicated. Source: `output/cgga_validation.json:cohorts`, ledger 2k. The IDH-wild-type primary analysis differs (ledger 2k); do not pool the populations. |

## Model
| ID | Limit | Status | Note |
|---|---|---|---|
| L8 | Adaptive results rest on assumed kill scale (E_MAX = rho x 1000) | **Checked: sensitivity 500-2000** | Simulation with assumed kill rates, ledger 2j: TTP split and rho > 0.02 split identical at 500/1000/2000; resistance advantage 57/54/51 of 61. Kill rates stay assumed. |
| L9 | Script 48 (3D) compared legacy and 10x-higher diffusivity regimes | **Checked, model sensitivity** | Alternate 0.13/0.013 regime changes selected high-rho results (0.678 to 0.346; 0.109 to 0.051); not proof of a literature-calibrated model or general robustness. See ledger 2h and `output/cohort_3d_summary.json`; regime comparisons are also summarized in `archive_old_gbm/vault/Script-48-3D-Extension.md`. |
| L11 | No clonal competition; resistance is a fixed fraction | Open | Model upgrade. |
| L12 | Anisotropy not matched to real tumour shape; alignment lower under aniso | Open | Later manifest contrasts show no measurable tract-anisotropy gain, ledger 2w/2y; the cause remains unproved. |

## Evaluation
| ID | Limit | Status | Note |
|---|---|---|---|
| L13 | Historical RL simulations with assumed kill rates: day-90 volume and TTP are distinct endpoints (ledger 2i) | **Finding** | |
| L14 | Pacing +2.1 d on real parameters | **Simulation finding** (see above) | |
| L15 | Probe-then-commit needs precise measurement | Open | No power analysis. |
| L16 | rho > 0.02 threshold: n=10, chosen after seeing data | **Checked** | 10/10 vs 0/51; permutation p 5e-5; leave-one-out 98.4%. Still post hoc. LUMIERE has only 3-4 patients above 0.02, so it cannot test it. |
| L17 | Forecast target chosen after a volume probe | Open | Script 81 `--target wt`, about 12 h. |

## Negative results
| ID | Limit | Status | Note |
|---|---|---|---|
| L18 | D_f interpretation withdrawn | Historical finding | Replaced by elongation; original validity analysis/source key is not retained in this summary. |
| L19 | MGMT does not predict TTP (n=130, p 0.60) | **Checked: underpowered** | 108 events; smallest detectable HR 1.74; observed 1.11 needs 2,807 events. CGGA OS is inconclusive (HR 0.91; pooled cohort n=374, MGMT-available n=328). Sources: `output/limit_checks.json:L19_mgmt_power`, `output/cgga_validation.json:cohorts.pooled_all_gbm.tests.mgmt_meth`. |
| L20 | Resistance-driven adaptive gain fails on real growth inputs under assumed kill rates | Finding | Scripts 76, 82, 83. |
| L21 | Historical growth-classifier AUC 0.46 | Unverified here | Exact output/key and endpoint are missing; do not cite this summary as evidence. |

## Not tested
| ID | Limit | Status | Note |
|---|---|---|---|
| L22 | No prospective validation | Not possible | |
| L23 | 3D alternate-D sensitivity (not established literature calibration) | **Fixed** (see L9) | |
| L25 | Drug screen is in silico only | Open | MT-CO2 is mitochondrial: nuclear CRISPR screens (DepMap) may not cover it. Check script 31 first. |

## New limits found in this round
| ID | Limit | Source |
|---|---|---|
| L26 | CGGA age effect is much weaker than TCGA and not significant after Holm in each batch | `output/cgga_validation.json:cohorts`, ledger 2k; c-index interpretation unverified |
| L27 | LUMIERE tumour volumes are noisy: median no-change error 0.8-1.8 log units; gap shrinks above 500 mm3 | `output/lumiere_volume_forecast.json`, ledger 2l; forecast error does not isolate measurement noise |
| L29 | RHUH-GBM has no MGMT column and only 40 patients | Historical metadata report; exact table/dictionary version missing here, unverified |
| L30 | Script 90 multivariate c-index is odd (pooled 0.48) and was not checked | `output/cgga_validation.json:cohorts`, ledger 2k; c-index interpretation unverified |
| L32 | The zone-expression files used by the stromal and adaptive scripts (43-44) are very likely synthetic (see ledger 2m) | `output/real_cohort_{le,ct,it}.csv` |
| L31 | Burdenko-GBM-Progression is restricted; not obtained | Historical image-access status unverified without dated approval/inventory; open clinical CSV exists (ledger 2az) |
| L33 | The "Healthy" class is 3 donors with other diagnoses (neurocytoma, oligodendroglioma, meningioma), all Fresh; no Frozen cell is Healthy | script 95, ledger 2o |
| L34 | MU-Glioma primary population mixes WHO grades: 30 of 134 eligible patients are not GBM | `data/manifests/consort_mu.json:mu`, ledger 2p |
| L35 | Eligible MU patients progressed far more often than excluded ones (88.8% vs 47.8%): follow-up is informative | `output/cohort_tables.json`, ledger 2p |

## Original work list (superseded)
Original work list superseded: all five script-88 folds and E_MAX sweep are reported above; RHUH later entered H-2/H-4. Other old proposals are not run authorization. Current RANO work and limitations are in START_HERE.md and ledger 2br–2bx.
