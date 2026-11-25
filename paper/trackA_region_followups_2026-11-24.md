# Track A region follow-ups (2026-11-24)

## Ravi et al. 2022 (Dryad doi:10.5061/dryad.h70rxwdmj), read directly

README (Dryad page, version 1.2, 2025-03-06): Visium samples are named by patient and region code ("C = Cortex, T = tumor, TC = Tumor Core, TI = Tumor Infiltration"). Each sample folder has Space Ranger outputs, `spatial_enrichment.csv`, and raw H&E images. No per-spot pathologist annotation file is listed. Region is a whole-section label, so region and section are confounded. CC0. Not usable for a pathologist per-spot test as released; usable if a pathologist annotates the included H&E images.

## Why the Saadawy edge signature ranked backwards (spent data, descriptive)

Checked on one spent section (AT4-BRA5-FO-1_1; 967 core, 1,114 edge spots). Not a code error: both halves of the signature are higher in pathologist-labelled core spots (leading-edge genes mean z +0.38 core vs -0.33 edge; core genes +0.29 vs -0.25). The leading-edge genes (e.g., PTPRZ1 2.27 vs 0.63, BCAN 1.31 vs 0.35, FABP7 1.65 vs 0.27 log-CPM) are malignant glial-progenitor genes, so at spot level they track tumor-cell density, which is higher in cellular tumor. A signature of edge-resident malignant cells therefore marks core tissue at spot resolution. Single section; descriptive only.

Implication: published "edge" signatures built from malignant cells should not be applied to mixed-tissue spots without adjusting for tumor-cell content.
