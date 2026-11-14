# Track A Wu specimen crosswalk and multi-study route

Read only specimen metadata columns from [Wu Supplementary Table 1](https://ars.els-cdn.com/content/image/1-s2.0-S1534580725002606-mmc2.xlsx). TC denotes tumor core, PT peritumoral brain, and TE tumor edge. For this screen PT is the distal margin candidate. Keep TE distinct and confirm the surgical definitions before freezing the endpoint.

| Primary IDH-wildtype origin | TC specimens | PT specimens profiled for RNA | RNA assay |
|---|---|---|---|
| 6237 | 6237_C, 6237_D | 6237_A | Evercode |
| 6425 | 6425_E | 6425_A, 6425_F, 6425_G | Multiome |
| 6467 | 6467_D, 6467_E | 6467_A, 6467_B, 6467_C, 6467_F | Multiome |
| 6509 | 6509_C, 6509_D | 6509_A, 6509_E | Multiome |
| 6514 | 6514_B, 6514_C | 6514_D | Multiome |

These five pairs are an upper bound before malignant-cell and shared-state QC. Source for every code and count in this table is the supplement above, joined by explicit specimen prefix to the RNA records. 6419 has TC and PT but is recurrent. 6245 has PT and TE, with no TC. These two facts explain why the author-named donor list from the [follow-up](trackA_fresh_data_followup_2026-11-13.md) cannot serve as the validation manifest. The two GEO RNA objects return identical Content-Length values, which do not prove identity. Keep both assays in the same study for validation.

## Other independent-source routes

The [Harwood CosMx study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12771808/) describes eight primary patients with pathologist-verified infiltrated tissue and core. Its [Zenodo record](https://zenodo.org/records/16957940) is open under CC BY 4.0. The central directories of its four archives list molecule tables and imaging files. The inspection did not identify a standalone tabular donor-region-malignancy crosswalk. Ask for FOV-to-donor and pathology labels plus malignant-cell annotations. The targeted panel also requires a fixed intersection and signature coverage audit. Harwood's pseudobulk method and tissue-stratified comparisons belong in the prior-art inventory. Screen overlap with the earlier Harwood dataset before assigning either to development or test.

Petti and Pai are candidates for distinct test sites if the crosswalk, diagnosis, comparable distal-region definition, and annotation gates clear. These labels must establish the endpoint independently from expression. The matched-margin protocol requires a fresh external study. A pooled multi-study test would require a dated amendment before expression access, with fixed studies, study-stratified reporting, one vote per donor, and no fitting on any test site. Pooling donors merely to cross the floor is insufficient. See [the existing cohort audit](trackA_alternate_dataset_plan_20261101.md).

## Pass-rule feasibility

For an exact one-sided sign test with n discordant donors, unanimous wins give p = 2^(-n). With m required comparisons, the smallest raw p must be at most alpha/m to pass the first Holm step. Freeze the comparator list, then check an attainable p-value bound and expected discordance before committing the final validation size. The [prospective protocol](trackA_matched_margin_protocol_2026-11-13.md) sets the significance threshold. No threshold was relaxed.

No candidate RNA matrix, molecule table, or drug response values were opened. The [candidate scoring code](../src/trackA_margin_pairs.py) was checked with synthetic arrays for consistent and inconsistent reversals, expression ties, no surviving pairs, monotone transformations, and gene-axis reordering. This verifies mechanics only. Published comparators and end-to-end input adapters are still pending.
