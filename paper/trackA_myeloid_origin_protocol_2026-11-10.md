# Track A exploratory analysis: myeloid origin program across Core and Periphery (protocol, frozen 2026-11-10)

**Status: FROZEN before any expression value is read for this analysis.** Exploratory. The GSE117891 cohort was already used by A5, so this is not fresh or independent confirmation. A5's code, result, and failed 5-of-6 rule are untouched.

## Question and justification

Within GBMap-annotated Myeloid cells, does a predeclared program contrasting microglia-homeostatic genes with blood-derived macrophage genes separate Core from Periphery sampling points, donor by donor?

Justification: the A5 within-Myeloid check (`paper/trackA_review_followup_2026-11-10.md`) found the A5 candidate near chance inside Myeloid cells while logistic regression was not. Published single-cell work reported that blood-derived macrophages are enriched in the tumor core and microglia at the periphery (Darmanis et al. 2017, Cell Reports; also Müller et al. 2017, Genome Biology; Bowman et al. 2016, Cancer Cell). This tests that published biology as a fixed, untrained score in a different cohort. It is a known hypothesis, so the novelty is low; the value is a donor-level test with matched-null and contamination controls.

Not a drug-response measurement. Labels are author-supplied MRI-marked sampling points, not cell origin.

## Data

`load_dataset()` from the sealed runner `src/189_trackA_yu2020_lodo.py` (same hashed inputs, same QC: ≥1,000 UMIs and ≥500 genes). Primary donors from `paper/trackA_yu2020/donors.tsv` (six strict donors). Cells with `cell_type == "Myeloid"`. Eligible donor: at least 10 Myeloid cells in each class. Expected from saved metadata counts: GS1, GS2, GS3, GS5, GS13 (GS6 has 5 Periphery Myeloid cells). If fewer than 4 donors are eligible, stop and report "not run".

## Features (fixed now; no gene selected from this data)

Expression: log1p(counts × 10,000 / library size). Each gene is z-scored within each donor's Myeloid cells (label-free; genes with zero variance in a donor contribute 0).

- Microglia-homeostatic set MG: P2RY12, TMEM119, CX3CR1, SLC2A5, P2RY13, GPR34, SALL1, OLFML3.
- Blood-derived macrophage set BMDM: TGFBI, ITGA4, CD163, ANXA1, LYZ, S100A8, S100A9, FCN1.
- Neoplastic-marker control set NEO (ambient tumor RNA check): SOX2, OLIG1, PTPRZ1, BCAN, NES, SOX4, CKB, VIM.

Gene lists are literature-derived markers chosen before scoring; an exact match to any one paper's published list is not verified. Genes absent from the matrix are dropped and reported. If MG or BMDM keeps fewer than 4 genes, stop and report "not run".

Program score per cell: S = mean z(MG) − mean z(BMDM). Prediction for Core: −S (hypothesis: Periphery cells have higher S). NEO control score: mean z(NEO), predicting Core.

## Endpoint and comparisons

Per eligible donor: ROC AUC for Core (label 1) vs Periphery (label 0) within Myeloid cells.

1. Program AUC (−S).
2. Matched random null: 1,000 random gene-set pairs drawn from genes detected in that donor's Myeloid cells, each pair matched to MG and BMDM in size and in mean-expression decile; null AUC computed the same way. Seed 20261110.
3. NEO control AUC.
4. Trained baseline: L2 logistic regression (C = 1, `lbfgs`, max_iter 2,000) on the 2,000 most variable genes of the training donors' Myeloid cells, standardized on training cells, leave-one-donor-out across eligible donors. AUC on the held-out donor.

Equal-donor means are reported. No cell-level p-value is used as donor-level evidence.

## Rule (frozen)

- **Positive:** program AUC > 0.5 in at least 4 of the eligible donors (of 5 expected), AND program AUC > the donor's 95th percentile of the matched null in at least 4, AND program AUC > NEO AUC in at least 4.
- **Negative:** program AUC ≤ 0.5 in at least 3 donors.
- **Inconclusive:** anything else.
- The logistic baseline comparison is descriptive (wins counted, no rule), because the program is untrained and the baseline is trained.
- If NEO AUC is above 0.5 in at least 4 donors, report that ambient tumor RNA or composition may contribute, whatever the outcome.

## Output

`output/trackA_myeloid_origin_20261110/result.json` (summary metrics only; no expression values). Script `src/222_trackA_myeloid_origin.py`. Weakest points to report: small Myeloid counts per donor (e.g., GS5 23 Core / 11 Periphery), exposed cohort, mixed-cell sampling labels, normalized UMI data without ambient correction.
