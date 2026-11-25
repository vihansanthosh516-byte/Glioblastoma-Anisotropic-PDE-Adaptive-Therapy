# Track A region: one labelled section to other sites of the same tumor (LOCKED, 2026-11-24)

Locked before any histopathology label value is read (only per-section presence was checked). Code `src/trackA_region_transfer.py` (hash in the result). One run: `py -3.14 src/trackA_region_transfer.py`. The script refuses to overwrite `output/trackA_region_transfer/result.json`. No retuning after the run. Supersedes the draft `trackA_region_one_section_transfer_plan_DRAFT_2026-11-24.md`.

## Claim tested (within-patient only)

From one pathologist-labelled Visium section, a model trained on that section ranks core vs infiltrating/edge spots on the same patient's sections from **other sampling sites**, better than every comparator, in every eligible patient.

## Data

GBM-Space Visium per-section AnnData (E-MTAB-17164). Labels: "percentage of spot overlap with IvyGAP histopathology annotations". Verified in the GBM-Space paper Methods (De Jong et al., Nature 2026, s41586-026-11234-8): "A pathologist annotated digitized images of the H&E-stained tissue sections"; "Two more pathologists independently validated the histological annotations."

Sections: those with any non-zero label (`histopath_coverage.tsv`) whose file name matches `ATxx-...-{site}_{replicate}.h5ad`. AT15 uses a different naming scheme and is excluded by this rule. Expected eligible donors (two or more labelled sites): AT10, AT3, AT4, AT13. AT6 has one labelled site and is excluded.

## Spot labels (fixed)

Core = cellular tumor overlap >= 50%. Edge = infiltrating tumor + leading edge overlap >= 50% and core < 50%. Other spots are not scored. If all overlap values in a section are <= 1, they are treated as fractions and multiplied by 100. A section is used only if it has >= 20 core and >= 20 edge spots.

## Pairs

Each eligible section is a source. Targets are eligible sections of the same donor from a different site. Donor score = mean AUROC over that donor's pairs.

## Expression features

log1p(CPM x 1e4) per spot; genes detected in >= 1% of spots of the section; 2,000 most variable genes in the source's scored spots. Target genes missing from the target are set to 0.

## Primary candidate

Within-section rank normalization (per-gene quantile across the section's scored spots, applied to source and target separately), then source-standardization, PCA (30 components, fit on source), balanced L2 logistic regression (C = 1).

## Comparators

1. Ablation: the same model without rank normalization.
2. kNN label transfer: 15 nearest source spots in source PCA space (alignment-free mapping).
3. Saadawy 2026 leading-edge vs cellular-core markers (untrained, `paper/trackA_artzi_signatures/c26_saadawy_le_vs_cellular_core.tsv`).
4. Bastola edge/core markers (untrained, `c17_bastola_edge_core_markers.tsv`).
5. GBM-Space spatial niche abundance, best single niche and sign chosen on the source.
6. GBM-Space cell-state abundances, logistic regression trained on the source.
7. PASTE2 alignment transfer: **not run** (package installation was not permitted in this environment).

## Pass rule

PASS only if at least 4 donors are eligible and, in **every** eligible donor, the candidate's donor-mean AUROC is higher than every comparator 1-6. Otherwise FAIL. With 4 donors no sign test can reach p < 0.05 (best one-sided p = 0.0625); this rule is a consistency rule, not a significance test.

## Secondary (descriptive, not part of the verdict)

Split-conformal sets (alpha 0.10) for the candidate, calibrated (a) on the source site's other replicate section and (b) in-sample on the source; report target-site coverage and singleton rate.

## Code checks (generated data only)

Planted signal with section-specific per-gene scale shifts: candidate beat all comparators in 4/4 donors (PASS). Random labels: FAIL. Loader tested on one real section for shapes only (no labels read).

## Weak points

- Four donors; within-patient claim only. Not generalization to new patients.
- PASTE2, the main published alignment comparator, was not run.
- The rank normalization is a known technique; novelty rests on the task framing (one section to distant sites of the same tumor, scored against pathologist labels), which a prior-art search did not find but did not exhaustively rule out (`trackA_few_patient_designs_prior_art_2026-11-24.md`).
- Ivy GAP categories may miss some structure (Greenwald 2024). Labels are spot-overlap percentages, so mixed spots are excluded at the 50% thresholds.
