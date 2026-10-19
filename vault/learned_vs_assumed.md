# Learned vs assumed ledger (masterplan §40, §96)

> **HISTORICAL BACKGROUND (corrected 2026-10-19):** Dated table of the older GBM platform (created 2026-10-04), not current status. The "fill the blanks" instruction and Phase 3-6 references belong to the old Track plan and are obsolete. Current goal: probabilistic RANO (`START_HERE.md`). Status words are informal labels without decision rules; "TBD" cells were never filled and those provenance claims are unverified.

Status words: Observed (in the dataset), Derived (computed from observed), Inferred (fitted), Assumed (chosen), Simulated.
Created 2026-10-04. (Old instruction "fill the blanks as Phase 4-5 runs" is obsolete.)

| Quantity | Status | How | Prior / value | External check? | Source |
|---|---|---|---|---|---|
| MRI tumour mask (MU) | Observed/derived | expert-refined segmentation | – | LUMIERE masks pending | `data/` |
| LUMIERE volumes | Derived | pyradiomics, two segmenters | – | – | script 91 |
| Scan dates | Observed | relative days | – | – | manifest TBD |
| Treatment dates | Observed where present | clinical table | missing rate TBD | – | manifest TBD |
| Growth rate rho | Inferred | inverse fit | 124/154 at lower bound (L5) | LUMIERE rates 0.0042 / 0.0007 /day | `51_*`, script 91 |
| Diffusion D | Inferred (weak) / Assumed | Swanson-regime fixed 0.13 mm²/d parallel, 0.013 perpendicular in re-runs (legacy values 0.013/0.0013; ledger 2e) | not identifiable (L5) | none | script 48, 78, 80 |
| Carrying capacity K | Assumed | fixed | TBD | none | PDE code |
| Anisotropy ratio | Assumed | 10× parallel vs perpendicular (synthetic) | – | DTI adds nothing (−0.001; development comparison only, external check not run) | Track B |
| DTI tensor | Derived | UCSF-PDGM tensors | – | – | `UCSF-PDGM-Tensor-Study` |
| Drug concentration | Assumed | model | – | none | scripts 44, 75 |
| Kill rate k | Assumed (real_test uses same k for 21 patients) | E_MAX/rho ratio 1000 in the baseline, ratios 500, 1000, 2000 swept; real_test kill medians 0.05/0.08/0.13 shared across 21 patients | – | none | script 44, 92, L14 |
| Resistant fraction f_r | Assumed | fixed fraction | – | none | script 44 |
| Resistance transition | Assumed | fixed, no clonal competition (L11) | – | none | script 44, 83 |
| Inflammation score | Assumed (=1.0 for all 61) | zone files likely synthetic (L32; unconfirmed) | – | none | scripts 43–44 |
| rho > 0.02/day cutoff | Derived post hoc | chosen after data | – | LUMIERE has 3–4 above (HD-GLIO-AUTO 4/58, DeepBraTumIA 3/64) | script 85, 87 |
| Treatment policy | Simulated | PPO, MPC, paced | – | – | `58_*`, 75 |
| Noise levels | Assumed | 0–20% | old-plan item; later segmentation-noise analyses exist (`START_HERE.md`) | – | – |
| Molecular programs | Inferred | TCGA/atlas | gene sets never named (TBD) | CGGA survival validation only, not a test of expression programs | scripts 36–37, 90 |
| Molecular → rho link | Hypothesised | none yet | – | – | Old-plan Phase 6 (as of 2026-10-04) |

Rule: if a row says Assumed, the paper must not call the result "learned" or "predicted".
