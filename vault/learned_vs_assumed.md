# Learned vs assumed ledger (masterplan §40, §96)

Status words: Observed (in the dataset), Derived (computed from observed), Inferred (fitted), Assumed (chosen), Simulated.
Created 2026-10-04. Fill the blanks as Phase 4–5 runs.

| Quantity | Status | How | Prior / value | External check? | Source |
|---|---|---|---|---|---|
| MRI tumour mask (MU) | Observed/derived | expert-refined segmentation | – | LUMIERE masks pending | `data/` |
| LUMIERE volumes | Derived | pyradiomics, two segmenters | – | – | script 91 |
| Scan dates | Observed | relative days | – | – | manifest TBD |
| Treatment dates | Observed where present | clinical table | missing rate TBD | – | manifest TBD |
| Growth rate rho | Inferred | inverse fit | 124/154 at lower bound (L5) | LUMIERE rates 0.0042 / 0.0007 /day | `51_*`, script 91 |
| Diffusion D | Inferred (weak) / Assumed | Swanson-regime fixed 0.13 mm²/d in re-runs | not identifiable (L5) | none | script 48, 78, 80 |
| Carrying capacity K | Assumed | fixed | TBD | none | PDE code |
| Anisotropy ratio | Assumed | 10× parallel vs perpendicular (synthetic) | – | DTI adds nothing (−0.001) | Track B |
| DTI tensor | Derived | UCSF-PDGM tensors | – | – | `UCSF-PDGM-Tensor-Study` |
| Drug concentration | Assumed | model | – | none | scripts 44, 75 |
| Kill rate k | Assumed (real_test uses same k for 21 patients) | E_MAX = rho×1000; range 500–2000 tested | – | none | script 44, 92, L14 |
| Resistant fraction f_r | Assumed | fixed fraction | – | none | script 44 |
| Resistance transition | Assumed | fixed, no clonal competition (L11) | – | none | script 44, 83 |
| Inflammation score | Assumed (=1.0 for all 61) | synthetic zone files likely (L32) | – | none | scripts 43–44 |
| rho > 0.02/day cutoff | Derived post hoc | chosen after data | – | LUMIERE has 3–4 above | script 85, 87 |
| Treatment policy | Simulated | PPO, MPC, paced | – | – | `58_*`, 75 |
| Noise levels | Assumed | 0–20% | to be justified from measured segmentation error (Phase 3.3) | – | – |
| Molecular programs | Inferred | TCGA/atlas | pre-declared gene sets TBD | CGGA partly | scripts 36–37, 90 |
| Molecular → rho link | Hypothesised | none yet | – | – | Phase 6 |

Rule: if a row says Assumed, the paper must not call the result "learned" or "predicted".
