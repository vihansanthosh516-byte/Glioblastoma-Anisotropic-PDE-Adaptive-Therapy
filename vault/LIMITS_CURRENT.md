# Limits: current status (2026-10-04, evening)

Single list for the Limitations section. The ideas for fixes are in [[LIMITS_FIX_OPTIONS]]. Numbers are in [[PAPER_FINDINGS_LEDGER]].
Status words: **Fixed** = rerun and the limit is gone. **Checked** = tested, the limit is now measured but still there. **Open** = nothing done. **Finding** = it is a result, not a bug.

## Biggest (change what the paper can claim)
| ID | Limit | Status | What we know now |
|---|---|---|---|
| L0 | Track A uses a random cell-level split, with label-derived graph edges and test-set epoch choice | **Checked; script 88 running** | Patient-ID lookup scores 82.8%, above C-GAT 78.7%. Patient-level LR 60.3%, RF 62.0%, scVI-LR 66.2% (script 86). Leak-free C-GAT: fold 0 only (clean PCA graph 77.1%, cVAE without label edges 75.7%, published edges 86.4%). Wait for all 5 folds. |
| L24 | The Track A to Track B link | **Checked: none** | Inflammation score is 1.0 for all 61 patients (IDs do not match). Track A stands alone. |
| L14/L28 | `real_test` is not a real-drug test | **Checked** | Real growth rates, but the same assumed kill rates for all 21 patients. Median kill/rho 77.7 vs 2.7-12.5 in synthetic sets (script 87). This explains the +2.1 d gain. |
| L10 | Script 47 dual-agent "rescue" | **Fixed (explained)** | MTD plus the same second drug: 358.9 d, AUC 1152. Dual-adaptive: 305.6 d, AUC 1015 (script 89). The second drug does the work. |

## Data
| ID | Limit | Status | Note |
|---|---|---|---|
| L1 | Forecast Dice only 0.26, median gain +0.005 | Open | LUMIERE masks (30 GB zip) not downloaded. |
| L2 | Forecast loses on shrinking tumours (53% shrank or stayed) | **Finding, replicated** | LUMIERE volumes: growth-only model loses to no-change (n=40 p 0.0003; n=50 p 0.004; script 91). |
| L3 | UCSF-PDGM: grade 2-3, one scan | Open | RHUH-GBM downloaded (n=40, 3 time points, segmentations), not analysed. |
| L4 | Only 21 real patients in Track C | Open | LUMIERE gives 58-64 growth rates (script 91) but they are not run through Track C. |
| L5 | D not identifiable; 124/154 growth rates at the lower bound | Open | |
| L6 | 15,000-cell subsample, one seed, no CI | Partly | Scripts 86/88 add 5-fold mean and SD. |
| L7 | TCGA n=150 with expression; age is known | **Checked** | CGGA GBM (n=374): age HR 1.009 (1.001-1.018), Holm p 0.08, vs TCGA 1.031. Sex, MGMT null (script 90). Gene models not replicated. |

## Model
| ID | Limit | Status | Note |
|---|---|---|---|
| L8 | Adaptive results rest on assumed kill scale (E_MAX = rho x 1000) | Open | |
| L9 | Script 48 (3D) used D 10x below literature | **Fixed** | Re-run at 0.13/0.013: result holds, sparing at high rho lower (0.678 to 0.346; 0.109 to 0.051). |
| L11 | No clonal competition; resistance is a fixed fraction | Open | Model upgrade. |
| L12 | Anisotropy not matched to real tumour shape; alignment lower under aniso | Open | Cause of the alignment result not checked. |

## Evaluation
| ID | Limit | Status | Note |
|---|---|---|---|
| L13 | RL wins only on day-90 volume, never on TTP | **Finding** | |
| L14 | Pacing +2.1 d on real parameters | **Explained** (see above) | |
| L15 | Probe-then-commit needs precise measurement | Open | No power analysis. |
| L16 | rho > 0.02 threshold: n=10, chosen after seeing data | **Checked** | 10/10 vs 0/51; permutation p 5e-5; leave-one-out 98.4%. Still post hoc. LUMIERE has only 3-4 patients above 0.02, so it cannot test it. |
| L17 | Forecast target chosen after a volume probe | Open | Script 81 `--target wt`, about 12 h. |

## Negative results
| ID | Limit | Status | Note |
|---|---|---|---|
| L18 | D_f invalid | Finding | Replaced by elongation. |
| L19 | MGMT does not predict TTP (n=130, p 0.60) | **Checked: underpowered** | 108 events; smallest detectable HR 1.74; observed 1.11 needs 2,807 events. CGGA OS also null (HR 0.91, n=374). |
| L20 | Resistance-driven adaptive gain fails on real data | Finding | Scripts 76, 82, 83. |
| L21 | Growth classifier failed (AUC 0.46) | Open | |

## Not tested
| ID | Limit | Status | Note |
|---|---|---|---|
| L22 | No prospective validation | Not possible | |
| L23 | 3D at literature D | **Fixed** (see L9) | |
| L25 | Drug screen is in silico only | Open | MT-CO2 is mitochondrial: nuclear CRISPR screens (DepMap) may not cover it. Check script 31 first. |

## New limits found in this round
| ID | Limit | Source |
|---|---|---|
| L26 | CGGA age effect is much weaker than TCGA and not significant after Holm in each batch | script 90 |
| L27 | LUMIERE tumour volumes are noisy: median no-change error 0.8-1.8 log units; gap shrinks above 500 mm3 | script 91 |
| L29 | RHUH-GBM has no MGMT column and only 40 patients | `data/external/rhuh/` |
| L30 | Script 90 multivariate c-index is odd (pooled 0.48) and was not checked | script 90 |
| L32 | The zone-expression files used by the stromal and adaptive scripts (43-44) are very likely synthetic (see ledger 2m) | `output/real_cohort_{le,ct,it}.csv` |
| L31 | Burdenko-GBM-Progression is restricted; not obtained | TCIA |

## Still to run
Script 88 (4 folds left); RHUH-GBM analysis; script 81 whole-tumour target; E_MAX sweep (L8); power analysis (L15); stronger forecast baselines.
