# Track A sealed protocol 3: TAM activation, tumour core vs 5-ALA invasive region (frozen 2026-10-30)

New, independent cohort: E-MTAB-8743 (10 IDH-wildtype GBM donors). The project has not used it before (checked
2026-10-30, `paper/open_data_readiness_addendum_2026-10-30.md`). Expression values were not read before this seal;
only the matrix header and sample annotations were checked against the SDRF.

## Why this list (fixed before this cohort)
Script 157 (UCSC atlas, 7 donors, exploratory) found core microglia higher in SPP1, APOE, CD163, CD74, HLA-DRA, C1QA
(`START_HERE.md`, script 157 line). The skeptic noted the list came from the same atlas (circular there). This test
asks whether it holds in an independent cohort whose region labels do not come from expression.

## Data (user-approved download, 2026-10-30)
- `data/external/E-MTAB-8743_author_matrix/Final_count_exon.csv` (2,580,293 bytes, SHA-256 483a759d...): authors'
  exon count matrix, 13,505 genes x 50 samples, from GitHub TheSYSTEMScellSIGNALLINGlab/5ALA-Infiltrative-GBM.
- `Annotation_col.csv` (SHA-256 506844519b...): sample labels; all 50 names match the official SDRF.
- Weak point: how the matrix was made (aligner, reference, gene filter) is not documented. Licence: MIT on the repo
  (code licence); no separate data licence found. Labelled exploratory-provenance.

## Samples used
Per donor: label "Core" (SDRF tumour core) and label "Inv" (SDRF invasive region, 5-ALA-guided). 10 donors, 20 samples.
Rim, Pos and Neg samples are not used.

## Preprocessing (fixed now)
log2(CPM + 1) from raw counts over the 20 used samples (library size = column sum of all 13,505 genes);
genes kept if non-zero variance across the 20 samples.

## Method (same as scripts 173/174 Track A)
zone_diffs from `src/173_trackA_senmayo_edge_core.py` with edge = "Inv": per-gene z-score across samples; per-donor
delta = Inv minus Core; E = mean delta of the list minus the median of 1,000 expression-decile-matched random sets
(seed 20261030). Exact two-sided sign-flip over 2^10 = 1,024 patterns on cohort-RMS-scaled E.
Adjusted test: E_adj = residual of no-intercept least squares of E on [H, N, M] (H = hypoxia E, N = neuron E,
M = myeloid-abundance E), then RMS scaling and sign-flip.

## Prediction and success (all required; alpha 0.05, two-sided)
Predicted direction: core > invasive (E < 0).
1. Unadjusted p < 0.05 and E mean < 0.
2. Adjusted (hypoxia + neuron + myeloid abundance) p < 0.05 and same sign.
If 1 passes and 2 fails: report "explained by hypoxia / normal-brain / myeloid content", not activation biology.
If 1 fails: report negative. No other list will be tested as a rescue on this cohort without a new, declared protocol.

## Lists (frozen files in this folder; hashes in SHA256SUMS)
- tam_activation_genes.txt: 6 genes from script 157 (above).
- hallmark_hypoxia_genes.txt, neuron_genes.txt: same files as protocol 2 (`paper/trackA_oxphos/`).
- myeloid_genes.txt: AIF1, CSF1R, PTPRC, ITGAM, CD14, TYROBP, FCER1G, LAPTM5. Pan-myeloid abundance markers chosen
  by us (ASSUMED list, not from one paper); covariate for how many myeloid cells a sample has.

## Known risks
- 10 donors: sign-flip p floor about 0.002; 3 covariates use 3 of 10 degrees of freedom; low power for the adjusted test.
- Core vs invasive differ in tumour, hypoxia, brain and immune content; covariates are partial controls only.
- Prior art: TAM differences between GBM core and periphery are published (e.g. macrophage-rich core vs microglia-rich
  periphery, Darmanis et al. 2017 Cell Rep). A positive result is a replication of a known pattern, not new biology.
  The authors of this dataset already compared core and invasive regions.
- Bulk tissue, not single cells: "activation" cannot be separated from cell-type mix beyond the covariates.

## Amendment 1 (2026-10-30, before any score was computed)
The first sealed run stopped with an error before any statistic was produced: the matrix rows are Entrez gene IDs,
not gene symbols, so no list gene was found. No expression value was summarized or shown. Fix: list symbols are
mapped to Entrez IDs with `symbol_to_entrez.tsv` (MyGene.info v3 query, human, exact symbol match; 283 of 283 mapped).
Lists, method, seed and success rule are unchanged.
