# Track A GBM-Space state and TMZ-response-associated score

**Status:** Protocol frozen for one exploratory analysis before reading or scoring the expression matrix. This analysis is separate from A5 and will not change the A5 model, endpoint, result, or figure.

## Question

Within the 12 primary IDH-wildtype GBM-Space tumors, does a fixed gene-expression score based on published in-vitro temozolomide (TMZ) viability associations differ between author-annotated malignant developmental-like states and author-annotated malignant gliosis/hypoxia states?

This is an exploratory association question. It does not test actual treatment response, predict a patient's response, or establish drug resistance. The directional gene panel was derived from in-vitro viability associations, not validated as a complete response predictor. A higher score means only that genes associated with higher in-vitro viability rank higher in that pseudobulk than genes associated with lower in-vitro viability.

## Prior-art boundary

The GBM-Space paper already establishes the developmental-to-gliosis/hypoxia malignant-state trajectory and its spatial organization. scFOCAL already integrates single-cell GBM states with drug-perturbation signatures and reports in-vitro, ex-vivo, and in-vivo validation. This project will not claim a new drug-response framework, a new state trajectory, or superiority to scFOCAL. The narrow exploratory contribution will be testing a separately published, fixed ex-vivo TMZ-viability-associated panel against the GBM-Space authors' state labels at the tumor level.

## Dataset and cohort

The analysis will use the public processed single-nucleus AnnData file `GBM_space_snRNA.h5ad` linked by the official GBM-Space portal. The paper reports 12 primary IDH-wildtype tumors, 57 sampled sites, and 1,025,329 nuclei. The local metadata-only audit found the same 12 donor IDs and 657,275 cells labeled `Malignant`; these counts and the donor-by-group table are recorded in `output/trackA_gbm_space_tmz_state/metadata_preflight.json`. The donor, site, `cell_status`, and `annotation_coarse` fields will be used; no nucleus will be treated as an independent patient.

The analysis will use the exact donor set frozen in the metadata audit: AT3, AT4, AT5, AT6, AT7, AT9, AT10, AT11, AT12, AT13, AT14, and AT15. Every donor must meet the cell-count rule below in both primary groups; no donor will be dropped to improve a result.

## Frozen state groups

Only cells with `cell_status == Malignant` and a nonmissing `donor_id` will be eligible.

- **Developmental-like:** `AC-progenitor-like`, `OPC-like`, `NPC-neuronal-like`, `OPC-NPC-like`, and `OPC-neuronal-like`.
- **Gliosis/hypoxia:** `AC-gliosis-like`, `Gliosis-like`, and `Hypoxic`.
- **Excluded from the primary contrast:** proliferative cells, TME cells, undefined labels, and any other state label.

The labels and grouping follow the source paper's annotation hierarchy. This is a comparison of author-defined state groups, not a re-annotation, spatial-region analysis, or causal state transition test. All sites will be pooled within each donor and state group. No site-level replicate will be counted as a donor.

Each donor must contribute at least 200 eligible malignant nuclei to each primary group. This threshold is fixed before expression scoring. The metadata-only counts exceed this threshold in all 12 donors. If a donor fails after the input is reopened, the complete primary analysis will stop rather than change the cohort.

## Fixed response-associated gene score

The source panel is the 121-gene `Correlation_cell_Overlap` panel from Ntafoulis et al. (2023), with directions defined by its reported correlation with in-vitro cell viability at 100 μM TMZ. The A5 protocol already maps 119 of these genes exactly to its count matrix. For a directly comparable gene set, this analysis will use the exact intersection of those 119 genes with GBM-Space `raw.var['SYMBOL']`.

Metadata-only mapping found 115 exact shared genes: 47 positive-viability (resistance-associated) directions and 68 negative-viability (sensitivity-associated) directions. No aliases, orthologs, substitutions, imputation, or outcome-based gene selection will be used. Six source genes will therefore be omitted: `AF131215.2`, `PPIAP22`, `CROCCP3`, `HERC2P2`, `HKR1`, and `LINCR-0001`. If the matrix's exact mapping or these counts differ at run time, the analysis will stop.

For each donor and state group, raw counts will be summed across eligible nuclei to form one pseudobulk. Pseudobulk counts will be converted to CPM, transformed as `log2(CPM + 1)`, and ranked across all measured features using average ranks for ties. Each feature rank will be converted to a percentile as `(rank - 1) / (G - 1)`, where `G` is the number of measured features. The score will be:

`mean(percentile rank of 47 positive-viability genes) - mean(percentile rank of 68 negative-viability genes)`.

The source matrix's `raw.X` will be used only if the following validation passes during a single sequential pass: the matrix has the frozen shape, every stored value is finite, nonnegative, and integer-valued within floating-point tolerance, and every row sum matches `obs.total_counts`. Failure stops scoring without writing a result. Expression-matching features will be estimated from a preselected simple random sample of 25,000 malignant nuclei, selected from metadata alone with seed `20261108`; the state-group labels will not be used to select matching cells.

## Primary comparison and controls

For each donor, calculate `gliosis/hypoxia score - developmental-like score`. The primary effect estimate will be the equal-donor mean of the 12 paired differences. Report every donor's two scores and paired difference, the mean, median, range, sign counts, and a 95% donor-level t interval. Also report an exact two-sided sign-flip p-value over all 4,096 donor-level sign assignments as an exploratory check. It will not be used as a pass/fail rule or as proof of a treatment effect.

Generate 1,000 expression-matched control panels of 115 genes. Match each source gene to a unique background gene among its 20 nearest available genes using mean per-cell `log2(CPM + 1)` and detection fraction, estimated on the same 25,000 malignant nuclei without using the two state-group labels. Standardize both features across measured genes; break distance ties by Ensembl ID; assign each control gene the source gene's direction; use NumPy `PCG64DXSM` seed `20261108`. The background will exclude every measured gene whose symbol appears anywhere in the original 121-gene source panel, including a source gene omitted from the 115-gene score. Each control panel will contain 47 positive and 68 negative directions. Report the primary donor-mean difference beside the median and 2.5th–97.5th percentiles of the control-panel donor-mean differences. These controls are descriptive and are not a second hypothesis test.

The primary contrast is the only inferential comparison. State subclasses and proliferative cells will not be tested separately in this run. No settings or gene definitions will be selected using the observed score.

## Stop rules and claim limits

Stop without a result if the local pre-score file fingerprint, matrix shape, donor set, required metadata fields, state labels, 200-nucleus donor-group threshold, 115-gene exact mapping, raw-count checks, or code/protocol seal fails. The committed fingerprint will record the local file size and SHA-256 before any expression scoring; it will establish which local bytes were scored, not independently authenticate them against an upstream publisher checksum. Do not silently remove donors, change groups, change the gene panel, or substitute another matrix slot after inspecting scores.

The result will be described as an exploratory association between a published in-vitro TMZ-viability-associated expression score and author-defined malignant-state groups in one 12-tumor atlas. Both the score and the state labels come from expression data, and the atlas has no matched TMZ response outcomes for these nuclei. The analysis will not establish that cells or patients are TMZ-resistant or TMZ-sensitive; predict clinical response; support a treatment recommendation; demonstrate causation; validate a full 121-gene signature; or claim an original drug-response method. The data will not be joined to Duke MRI patients or used as a patient-specific Track C treatment parameter.

The article is CC BY 4.0 and the official portal distributes the processed atlas through CELLxGENE. The local matrix does not itself carry a license field. The source matrix and cell-level derivatives will remain local; any public summary will include attribution to the GBM-Space source and will be limited to aggregate donor-level results. The exact dataset reuse terms will be checked again before external redistribution.

## Provenance and output

Before scoring, the local input fingerprint will be recorded in `paper/trackA_gbm_space_tmz_state_input_fingerprint_2026-11-08.json`, and the code, protocol, input fingerprint, and source scoring files will be sealed in `paper/trackA_gbm_space_tmz_state_SHA256SUMS`. The runner will verify the input filename, byte count, and SHA-256 against that frozen fingerprint before opening the AnnData file. The run will record SHA-256 hashes for the input matrix, source panel, this protocol, code, and the relevant A5 score implementation; the 12 donor IDs; group cell counts; feature mapping; deterministic seeds; package versions; and runtime environment. It will save aggregate donor-level scores and the control summary to `output/trackA_gbm_space_tmz_state/result.json`, refuse to overwrite an existing output, and leave the source matrix local. No project figure will change in this analysis.

## Sources

- de Jong, G., Memi, F., Gracia, T., et al. (2026). A spatiotemporal cancer cell trajectory underlies glioblastoma heterogeneity. *Nature*. https://doi.org/10.1038/s41586-026-11234-8
- Suter, R. K., et al. (2026). Drug and single-cell gene expression integration identifies sensitive and resistant glioblastoma cell populations. *Nature Communications, 17*, 99. https://doi.org/10.1038/s41467-025-67783-5
- Ntafoulis, I., et al. (2023). Ex vivo drug sensitivity screening predicts response to temozolomide in glioblastoma patients and identifies candidate biomarkers. *British Journal of Cancer, 129*(8), 1327–1338. https://doi.org/10.1038/s41416-023-02402-y
