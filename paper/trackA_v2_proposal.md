# Track A v2: provisional proposal

Status: research only. No patient data downloaded or scored. No dataset has yet passed the main-test eligibility gate. The prior-art search and download inventory are incomplete; do not seal or run this proposal yet.

## Dataset candidates

1. **GBM-Space**, [E-MTAB-17183](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-17183), [atlas](https://www.gbmspace.org/). Atlas reports 12 GBM tumors with multiple sampling sites. Tissue histopathology is available, but whether independent core/edge labels and both regions exist for at least 10 distinct donors, and whether single-nucleus samples map to those labels, is unverified. Data terms and exact processed file names, URLs and sizes are unverified. Raw data include controlled-access EGA assets. Best candidate for further read-only checking, not yet eligible. Weakest point: donor count is not the paired eligible donor count.
2. **Regional GBM**, [EGAS50000000547](https://ega-archive.org/studies/EGAS50000000547), dataset EGAD50000000778. The archive describes 11 regionally sampled tissues, not 11 donors. Sampling covers core, edge and peritumoral region; independent label procedure and sample-to-donor mapping remain unverified. EGA data require access approval; exact file names and sizes require catalog review. Weakest point: eligible donor count is unknown.
3. **GSE162631**, [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162631): four donors with matched surgically sampled tumor core and peripheral tissue. GEO lists `GSE162631_raw_counts_matrix.zip.gz`, 806.1 MB, at https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE162631&format=file&file=GSE162631%5Fraw%5Fcounts%5Fmatrix%2Ezip%2Egz . GEO terms and original authors' reuse conditions need review. Weakest point: below the ten-donor minimum and endothelial-enriched; cannot be the main test.

Counts and file sizes above are repository claims, not local output results. They have not been independently audited. File inventory is incomplete and is not a download request. No data download is authorized by this document.

## Method proposal

Use [scVI](https://doi.org/10.1038/s41592-018-0229-2), a neural generative model of single-cell counts. The current [scvi-tools license](https://github.com/scverse/scvi-tools/blob/main/LICENSE) is BSD 3-Clause. Use CPU and verify dependency compatibility before implementation. scVI is not a graph network.

One modification: a donor-adversarial loss on the latent representation, trained only against donor IDs in the training folds. It may suppress donor-specific shortcuts; it may also suppress genuine biology. Test this, never assume success. A shared downstream region classifier uses the same donor-balanced training, inputs and prediction aggregation for both models. Set adversarial weight to zero for the original scVI ablation. Held-out donor identity is never supplied to training; fit normalization, feature choice, model and classifier on training donors only, including each inner fold. Avoid transductive fitting on test donors.

Primary metric: mean of within-donor region AUC for donors with both independent labels and sufficient cells in each region. Choose minimum cells and class balance before sealing. Compare original and improved model on identical nested donor-held-out folds via `src/donor_cv.py`. Baselines: training prevalence, cells-per-donor, cell-type composition, and regularized pseudobulk gene expression. A cells-per-donor-only predictor is constant within each donor, so its within-donor AUC is 0.5; also report donor-level region classification where that baseline is meaningful, without mixing estimands. Use paired donor bootstrap intervals on out-of-fold donor scores, with training variation acknowledged.

Bar 1 passes only if the lower bound of the prespecified paired donor bootstrap interval for improvement over original scVI exceeds zero, the improved model exceeds every valid baseline on the same donors and metric, and the ablation returns to the original result. A separate donor-level endpoint must beat the cells-per-donor baseline with a paired interval above zero. Specify endpoint and interval details before sealing. Repeated tuning, small donor counts, and overlapping training sets weaken uncertainty claims. No superiority is established.

## Candidate biology question

Among malignant cells from independently sampled core and edge, does the association between unfolded-protein-response activity and interferon response differ by region within donors, after adjustment for hypoxia and normal-brain contamination? This is a candidate, not a novelty claim. Define independent external gene sets, a paired donor estimand and covariates before data access. Avoid circular analysis: region labels must come from tissue sampling, not gene expression. Specify how model-derived cell states are used and preserve an independently computed program score as a sensitivity analysis. Require a prespecified association interaction with a donor-level uncertainty interval excluding zero in both unadjusted and confound-adjusted analyses; test donor composition by within-donor contrasts and leave-one-donor-out sensitivity. Malignant identity and contamination need independent validation. Do not report a positive biology claim if any gate fails.

Closest located papers, subject to further checking:

- [Pastore and De Rango, 2026](https://pubmed.ncbi.nlm.nih.gov/42214146/) study an interferon, stress and unfolded-protein response axis across GBM samples. Their spatial strata overlap this idea; the proposed paired region interaction is narrower. Whether already tested requires full-text review.
- [Pai et al., spatial multi-omics of the margin](https://pmc.ncbi.nlm.nih.gov/articles/PMC12603948/) compare matched core and margin and describe infiltrative programs. Proposed interaction differs, but their full results and supplements must be searched before claiming novelty.
- [Lim et al., 2026](https://pubmed.ncbi.nlm.nih.gov/41516347/) study necrosis effects on IRE1-XBP1 signaling in GBM cell lines. Proposed analysis concerns a donor-paired regional interaction in human tissue, with hypoxia adjustment.

Search status: PubMed-indexed papers located by web search. Direct PubMed, Europe PMC and Google Scholar query pages were inaccessible in this session. Consequently this is not the requested completed three-database prior-art check. Repeat all three searches and review full texts and supplements before building; reject this question if already answered.

## Gates before any download

Resolve candidate donor/sample/region mapping, independent label provenance, reuse terms, and list each required file with exact URL and verified size. Present that list for explicit user approval. Do not download until approval. Before scoring, write and hash the protocol, gene lists and settings; commit and push the seal, following `paper/trackA_senmayo/` and `src/173_trackA_senmayo_edge_core.py`. Record all results in `output/`, include negative results and the weakest point, and label exploratory analyses. No clinical benefit claim.
