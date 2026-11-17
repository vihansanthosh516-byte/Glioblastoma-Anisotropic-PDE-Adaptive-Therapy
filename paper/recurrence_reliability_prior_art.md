# Recurrence measurement reliability: closest prior art and novelty gate

Date: 2026-11-16. Phase A of `paper/CLAUDE_RESEARCH_PAPER_HANDOFF_2026-11-16.md`. Literature only; no expression values were opened for this note.

## Decision

**Gate passes, narrowly.** I found no study that asks how often a patient's primary-to-recurrence direction for one fixed malignant-cell program changes when only technical processing changes. The closest studies either (a) estimate patient-level recurrence trajectories without testing technical reliability (Alexopoulos and Karampelas 2026), (b) benchmark scoring methods at the cell or case-control level without paired patient direction (Ciernik et al. 2026; PathwayBench 2026; Noureen et al. 2022), or (c) show that composition can dominate mesenchymal scores (Kreis et al. 2024). The contribution is therefore a **measurement-reliability method and benchmark**, not new biology. It must not be described as a discovery of recurrence biology, a new recurrence program, or external replication of Wang et al.

Closest competitor: Alexopoulos and Karampelas (2026), for data and patient-level framing. Closest methodological competitor: PathwayBench (2026), for donor-level scoring robustness.

What is unverified: the full texts of Spitzer et al. and Nomura et al. (open in PMC, abstracts only read here), the full methods and supplement of Alexopoulos and Karampelas (paywalled; abstract and data statement read), and the full texts of PathwayBench, the 465-pipeline normalization study, and sctrial (abstracts only). A study in press or outside the indexes searched could overlap. Searches: Firecrawl research index (PubMed, bioRxiv, medRxiv, arXiv) with three focused queries, plus a web search for the closest paper, on 2026-11-16.

## Close studies

| Study | Evidence read | Question | Data / patients | Program and scoring | Patient-level uncertainty | Malignant annotation | Technical-reliability test? |
|---|---|---|---|---|---|---|---|
| Wang et al. 2022, Nat Cancer, doi 10.1038/s43018-022-00475-x | Full text (PMC9767870) | Recurrence biology under standard therapy | 86 snRNA specimens, 49 patients; 36 primary/first-recurrence pairs | PN-MES axis by multiple-factor analysis / PCA; Seurat clustering; AddModuleScore for immune states | Cohort-level tests | CONICSmat CNV calls (LRT < 0.05, BIC difference > 50, posterior > 0.5) | No |
| Spitzer et al. 2025, Nat Genet, doi 10.1038/s41588-025-02168-4 | Abstract (research index; open full text PMC12081298 not read for this note) | Longitudinal ecosystem trajectories (CARE) | 59 patients, matched primary/recurrent; snRNA + bulk DNA | Malignant states and composition change | Patient trajectories (abstract level); technical reliability not mentioned | Not verified here | Not in abstract |
| Nomura et al. 2025, Nat Genet, doi 10.1038/s41588-025-02167-5 | Abstract (research index; open full text PMC12081307 not read for this note) | Multilayer transcriptional architecture | 121 samples, 59 patients (same CARE material as Spitzer) | New program architecture | Not this question | Not verified here | Not in abstract; supports that constructs differ |
| Alexopoulos and Karampelas 2026, OMICS, doi 10.1177/15578100261499487 (online 2026-10-25) | Abstract, data statement, references (full text paywalled) | Patient-level recurrence routes; candidate EGFR/NFIA module | **9** Wang pairs; 59 CARE patients | Recurrent-minus-primary Neftel program vectors; 10-gene module selected in Wang, tested in CARE | Exact proportion CI, sign-flip test, CI on group difference | Wang "tumor-classified" cells; CARE authors' CNV table | **No** technical-procedure reliability reported in the abstract |
| Ciernik et al. 2026, Genome Res, doi 10.1101/gr.280926.125 (ANS) | Abstract | Stability and comparability of signature scores | 9 scRNA data sets | Seurat, Scanpy, UCell, JASMINE vs new ANS | Cell-level | Not central | Yes, but cell-level annotation, not paired patient direction |
| PathwayBench 2026, doi 10.21203/rs.3.rs-10271634/v1 | Abstract | Donor-level pseudobulk pathway scoring robustness | 8 case-control sets, 682 donors | ssGSEA, GSVA, z-score, AUCell, UCell | Donor-level, unpaired case-control | Not tumor | Yes, but unpaired groups and not malignant programs |
| Noureen et al. 2022, eLife, doi 10.7554/eLife.71994 | Abstract | Benchmark of signature scoring in cancer scRNA | Several cancer sets | Bulk vs single-cell scorers | Cell-level | Cancer vs normal | Yes, cell-level only |
| Kreis et al. 2024, Cancer Res Commun, doi 10.1158/2767-9764.CRC-23-0383 | Abstract | Do stromal cells drive EMT/MES scores? | Bulk, scRNA, pseudobulk | 11 EMT/MES signatures | Not paired | Cell-type composition | Composition, not technical procedure |
| Johri et al. 2025, Nat Commun, doi 10.1038/s41467-025-57377-6 | Abstract | Group-level signature enrichment with hierarchical structure | Clinical scRNA | Group biology estimation | Group robustness | Tumor compartment | Hierarchical inference, not paired direction |
| sctrial 2026, bioRxiv doi 10.64898/2026.04.02.716219 | Abstract | Participant-level longitudinal estimands | Simulation and trials | Gene-level DE | Participant-level | Not central | Inference method, not scoring reliability |
| Saxena 2026, Comput Biol Chem, doi 10.1016/j.compbiolchem.2026.109100 | Abstract | Effect of 465 pipelines on marker interpretation | Cancer scRNA | Biological Discordance Score | Not patient-paired | Not central | Pipeline-level, not patient direction |
| Varn et al. 2022, Cell, doi 10.1016/j.cell.2022.04.038; Hoogstrate et al. 2023, Cancer Cell, doi 10.1016/j.ccell.2023.02.019 | Abstracts | Bulk longitudinal glioma | GLASS cohorts | Bulk deconvolution | Patient pairs | Bulk | No; composition explains part of recurrence change |

## What is new and what is not

Not new: pairing Wang patients; a mesenchymal shift at recurrence; adding CARE; patient-as-unit inference; showing that scoring methods disagree at the cell level; showing composition affects mesenchymal scores.

Candidate new element: a patient-level **direction-reliability estimand** for paired recurrence data. For each patient it asks whether two fixed procedures for the same frozen program agree on the sign of the recurrent-minus-primary change, with an indeterminate zone set by a declared within-specimen cell-resampling model rather than by the recurrence results. It reports opposite-direction discordance, indeterminate coverage, and continuous change together. I found no prior paper with this estimand. This is a narrow methods contribution. Its value depends on whether the result is informative at the attainable sample size (Phase D).

## Consequence for the protocol

1. Frame the study as measurement reliability of patient-level conclusions in the Wang reanalysis cohort.
2. Cite Alexopoulos and Karampelas as the closest data competitor. Do not claim to be first to study patient-level recurrence trajectories in Wang.
3. Treat construct disagreement (Neftel MES vs Verhaak mesenchymal) as secondary and not as technical error (Nomura et al. 2025 supports that constructs differ).
4. Published finding known before scoring: Wang et al. report a cohort-level mesenchymal shift. This is literature knowledge, not project outcome exposure.
