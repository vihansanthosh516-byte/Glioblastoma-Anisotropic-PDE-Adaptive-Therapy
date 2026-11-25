# Track A region development: why confidence breaks across sites (2026-11-24)

DEVELOPMENT on spent data (AT10, AT3, AT4 were spent by the locked test). Code `src/trackA_region_dev_conformal.py`; output `output/trackA_region_dev_conformal/dev_result.json`. 34 source-target pairs, conformal level 90%.

| Calibration | Mean target coverage | Worst pair | Singleton rate |
|---|---:|---:|---:|
| Marginal (replicate section) | 0.590 | 0.140 | 0.824 |
| Class-conditional (Mondrian) | 0.479 | 0.049 | 0.727 |
| Label-shift weighted (EM prior, Podkopaev-Ramdas style) | 0.484 | 0.000 | 0.654 |

- The share of core spots differs a lot between sites (mean absolute shift 0.35).
- The standard label-shift prior estimate (Saerens EM) was further off than doing nothing (mean error 0.40 vs the 0.35 shift), so label-shift corrections made coverage worse.
- Interpretation: the shift is not only in how many core spots there are; the score distribution of each class also moves between sites (conditional shift). Ranking (AUROC 0.90-0.97) survives; probability calibration does not, and textbook fixes do not repair it.

Fresh labelled data search (helper agent, 2026-11-24): no other public GBM spatial dataset with pathologist per-spot labels and >= 2 annotated sections per patient was found. Sonpatki 2026 (Zenodo 17572905, 17 patients, 32 slides) has per-spot anatomical labels, but they are derived from cell-state scores with manual thresholds and label propagation, not a pathologist. Ravi 2022 (Dryad) annotator and per-spot format unverified. The only untouched pathologist-labelled sections in hand are GBM-Space AT15 (one donor, CytAssist v2, at least 4 labelled tissue blocks; 4 files still downloading).

Next method idea (not yet built): target-adaptive calibration that uses the unlabeled target section's own score distribution (e.g., a two-component mixture on target scores, or spatial smoothness) instead of source-side calibration. Develop on spent data; test once on AT15 (fresh, but one patient) and, as a separate check, on Sonpatki labels (not pathologist).
