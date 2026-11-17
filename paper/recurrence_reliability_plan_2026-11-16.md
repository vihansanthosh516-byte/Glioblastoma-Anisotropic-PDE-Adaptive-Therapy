# Recurrence measurement reliability: feasibility and draft plan

Date: 2026-11-16. Research planning only; this is not an adopted protocol or authorization to open expression values. Existing A/B/C gates and reserved cohorts remain unchanged.

## Recommendation

Develop a focused measurement-reliability study before proposing another recurrence predictor. Ask whether a patient's inferred primary-to-recurrence change remains stable under technically defensible ways of processing the same measurements. This extends Track A with longitudinal tissue evidence; it does not provide a patient join to Track B or treatment-response calibration for Track C.

The data-file listing gate now passes for the conservative 27-pair candidate subset. The scientific eligibility and novelty gates do not yet pass. The main remaining requirements are reliable malignant-cell annotations, library/batch interpretation, first-recurrence provenance, and a contribution distinct from recent recurrence studies.

## What file metadata establishes

The [official GEO file list](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174554/suppl/filelist.txt) contains a complete matrix/features/barcodes file triple for each of the 54 candidate specimens. There are 61 complete library triples, because seven specimens have multiple listed libraries. Those libraries are not extra patients. The listed candidate files total 713,265,211 compressed bytes. The enclosing `GSE174554_RAW.tar` is 1,655,029,760 bytes. These are archive-listing sizes, not a verified minimal independent download route. No matrix, feature list, or barcode contents were opened.

Multiple-library specimens: SF10099, SF12243, SF4449, SF4810, SF7307, SF10433, SF4209. Before analysis, establish whether these represent disjoint captures, resequencing, or another batch relationship. Do not concatenate or sum them blindly, and do not count libraries as biological replication.

See [processed coverage manifest](../output/paper_strengthening_metadata/recurrence_processed_coverage.json), [source file list](../output/paper_strengthening_metadata/recurrence_filelist.txt), and [specimen crosswalk](recurrence_specimen_crosswalk_2026-11-16.md). File presence does not establish integrity, raw counts, suitable cell quality, annotation validity, or usable malignant-cell support.

## Closest prior work changes the proposed contribution

- Wang et al. already reports paired recurrence biology, including a mesenchymal shift. A new analysis of that same cohort is reanalysis, not external replication. [Original study](https://www.nature.com/articles/s43018-022-00475-x).
- Spitzer et al. analyzes longitudinal ecosystems from 59 patients, with cell-state and microenvironment composition changes. It is close prior art for a generic recurrence/composition question. [Nature Genetics 2025](https://doi.org/10.1038/s41588-025-02168-4).
- Alexopoulos and Karampelas, published online October 25, 2026, explicitly reanalyzes patient-level recurrence trajectories in Wang and CARE. Its accessible abstract also discusses trajectory uncertainty and malignant-cell annotation. Therefore pairing patients, comparing recurrence programs, or adding CARE alone does not establish novelty. Full methods were restricted during this review; the exact overlap remains to be checked. [Primary publication and accessible abstract](https://doi.org/10.1177/15578100261499487).
- Nomura et al. expands the architecture of GBM cellular programs. Alternative classification schemes need not measure the same construct; their disagreement cannot automatically be labeled measurement error. [Nature Genetics 2025](https://doi.org/10.1038/s41588-025-02167-5).

The possible contribution is a reproducible assessment of how technical choices affect patient-level conclusions, with an explicitly separated comparison of genuinely different biological state definitions. This is a candidate contribution, not a confirmed novel method. Do not launch the study as a new biological discovery unless the closest-method review supports that claim.

## Draft research question and estimands

Proposed question: How often does a patient-level recurrence-state conclusion change when technical sampling, library handling, or scoring changes, while the biological target is held fixed?

Choose one published malignant-cell program before viewing expression outcomes. The primary analysis compares two prespecified scoring procedures for the same frozen genes and intended construct. Compare different published state definitions only as a secondary construct-sensitivity analysis; there is no assumption that unlike states should agree.

The unit of inference is the patient. Define a continuous recurrent-minus-primary score for every patient. A primary reliability endpoint could be the proportion of patients whose directional conclusion changes across the fixed technical procedures. Report the continuous differences and their uncertainty as well; a small change around zero is not strong biological reversal. Define an indeterminate zone before scoring, using external technical evidence or a declared measurement model. Do not tune that zone to obtain a favorable flip count.

Nested cell resampling can describe within-specimen sampling variation but cannot substitute for patient replication. A patient bootstrap or a declared paired analysis addresses between-patient variability. Technical scoring replicates, cells, regions, and libraries are not independent patients.

## Comparators and controls

1. A simple frozen gene-set score on patient/specimen pseudobulk is the reference; compare a more complex scorer only if it targets the same construct.
2. Compare pooled-cell and equal-patient summaries descriptively to expose weighting effects. They estimate different quantities and are not competing evidence for the same population target.
3. Analyze author/CNV-supported malignant cells separately from whole tissue. Whole-tissue versus malignant-cell differences can reflect composition; they do not by themselves establish a malignant-state change.
4. Perform prespecified matched-depth cell subsampling and batch-aware sensitivity checks only after library relationships are established. Derive sampling-depth requirements before inspecting recurrence-score differences.
5. Use a finite, frozen list of procedures. Do not search many pipelines and retain those that agree. Report every planned procedure, exclusions, and failure.

## Precision and failure rules

There are only 27 metadata-matched candidate patients before biological QC. Even with zero observed failures among 27 independent patients, the exact one-sided 95% binomial upper limit on a failure probability is `1 - 0.05**(1/27)`, approximately 10.5%. Therefore this cohort alone cannot establish a failure rate below 10% at that confidence level, even in the best observed case. It may be smaller after QC. This calculation assumes independent patients and a fixed, consistently defined binary endpoint; technical resampling does not increase its denominator.

Do not choose a reliability threshold merely because 27 patients can pass it. Select a scientifically meaningful tolerance and calculate the required independent patient count before locking the protocol. If the attainable precision is insufficient, label the analysis exploratory or seek eligible independent data. No positive result is promised.

Stop before scoring if malignant-cell labels cannot be linked to specimens without circular target construction; candidate matrices lack usable support; batch relationships remain unresolved; or the primary endpoint remains indistinguishable from already-published work. Stop the confirmatory claim if eligibility, precision, or external replication fails. Preserve the negative result rather than changing endpoints after inspection.

## External validation lead, not a verified cohort

The recent recurrence paper points to CARE accessions GSE274546/GSE274548 and a public malignant-annotation repository. Official series text currently lists 121 GSMs for GSE274546 (10X) and one GSM for GSE274548 (Smart-seq2). These are accession/assay counts, not verified paired-patient counts. The 59-patient figure is a published study count. Confirm paired identifiers, modality, origin overlap, annotation linkage, rights, and prior project exposure before reserving any source as external validation. Two accessions may describe related material and are not automatically independent cohorts.

[Official GSE274546 metadata](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546&targ=self&form=text&view=brief), [GSE274548 metadata](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274548&targ=self&form=text&view=brief), and [local metadata extraction](../output/paper_strengthening_metadata/care_series_metadata.json).

## Execution order

First finish metadata, annotation, batch, access, and closest-method checks. Then write and commit a complete protocol, including the fixed gene program, scorers, indeterminate zone, QC, eligibility list, precision target, multiplicity, and failure rule. Only then inspect authorized expression data and run the declared analysis. Reserve external data before using it for tuning. Package the study around one reliability claim; A/B/C integration remains separate and simulation-only where applicable.
