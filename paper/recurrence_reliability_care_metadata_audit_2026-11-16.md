# CARE metadata eligibility audit (Phase G)

Date: 2026-11-16. Metadata only. No expression matrix, gene list, outcome, cell-state result, or score opened. GEO brieftext records and two small GitHub metadata files were inspected. Curated record: `output/recurrence_reliability_care_metadata/care_metadata_summary.json`.

## Decision

**CARE is a conditional metadata lead; do not download expression matrices or count it as 59 independent eligible patients yet.** Official records and the source paper support 121 10x tumors from 59 paired patients, with primary and first-recurrence labels. These are nominal pairs, not patients passing this project’s cell-label, barcode-join, and malignant-support criteria. Current CARE eligible independent count: **0 confirmed**.

GSE274546 brieftext contains 121 GSMs and 59 unique `Patient N` labels. GSM titles identify Timepoint 1/2; source labels say Primary/1st Recurrent. GEO characteristics show both Protocol 1 and Protocol 2. Spitzer et al. report 121 tumors/59 patients, 429,305 10x cells passing their QC, and seven collection institutions: MD Anderson, Duke, Tokyo University, Pitié-Salpêtrière, St. Michael’s Hospital, Seoul National University, and NORLUX. They classify malignant/nonmalignant cells using inferred CNAs. This gives a sound nominal paired design, but it does not establish CARE barcode-to-GSM-to-specimen malignant-label linkage or each specimen’s usable malignant-cell count. [GSE274546](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546), [Spitzer et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12081298/)

## Platform and patient independence

The companion paper reports that 54 of 121 10x tumors were also profiled using Smart-seq2. GSE274548 itself has one GSM, titled “GBM single-nucleus, Smart-seq2,” with a metadata field of 14,401 nuclei. The GEO record does not expose patient/timepoint identities for those nuclei. Treat Smart-seq2 as an overlapping assay subset, never 54 extra patients; it adds zero verifiable paired patients from current metadata. [GSE274548](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274548), [Nomura et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12081307/)

CARE source institutions listed in the Spitzer paper do not include UCSF, while Wang GSE174554 uses UCSF `SF...` specimen IDs. This source-level provenance supports CARE as a distinct cohort. It is not a patient identity crosswalk, so do not add 59 to the locked Wang cohort as independent n until the CARE-to-Wang identity question is closed with source-level linkage. The GitHub repository includes a separate 68-row `diaz_clinical_snrna_metadata_20240106.txt` with SF IDs and 34 pair labels; its SF IDs match the Wang comparison cohort, not CARE’s DU/SN GEO sample identifiers. It must not be mistaken for CARE’s 10x annotation or evidence of CARE/Wang overlap. [GBM-CARE-WT repository](https://github.com/dravishays/GBM-CARE-WT)

## Annotation and barcode gate

Paper methods say Cell Ranger v3.1.0 generated 10x count matrices; broad cell type labels were based on marker genes and inferred CNA classification. This describes the biological labeling method, not the needed public join implementation. The small GitHub linker `care_wt_synapse_dna_snrna_id_linker.txt` maps `P#T#` patient/timepoint labels to CARE DNA aliquot barcodes; it does not map a GEO GSM to 10x cell barcodes or malignant labels. The repository’s 10.9 MB `care_full_gene_count_metadata.RDS` was not opened. The separate `diaz_clinical_snrna_metadata_20240106.txt` is sample-level clinical metadata for SF comparison samples, not per-cell CARE annotations.

Before eligibility, verify a small metadata artifact/schema that joins each 10x barcode to GSM/specimen, patient, timepoint, and author CNA/malignant label. Then apply the locked protocol’s ≥50 malignant nuclei per specimen and technical QC gates. The paper’s 429,305 pooled cells do not establish per-pair support.

## Access and terms

Processed 10x data are listed on GEO as public; no individual access request is indicated for those processed files. Raw sequence reads are limited by consent and require controlled DUOS access, so they are outside this route. De-identified mutation/CNA tables are linked through Synapse; access and terms for exact files remain unverified. GEO public availability does not itself establish permission for the intended analysis. The paper’s CC BY-NC-ND article license does not license deposited data. **Access class: public processed-file route; terms: unverified at file level.** Keep full expression files unopened until exact terms and required annotation joins pass. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546), [Synapse CARE](https://www.synapse.org/care_glioblastoma), [DUOS](https://duos.broadinstitute.org/)

## Protocol compatibility and next gate

The locked Wang scoring reference estimates bins, control genes, and centering from all malignant nuclei in the final cohort. CARE must therefore use a prespecified transported Wang reference or a rebuilt CARE-specific reference; the latter changes the reference pool and is not an identical scorer deployment. Also, exact binomial inference assumes independent patient indicators; shared reference-pool estimation should be acknowledged and assessed before calling an external reliability estimate confirmatory. With 59 patients and zero discordances, the two-sided exact 95% upper bound is about 6.1%, still above the locked 5% tolerance.

Next allowed work: resolve processed-file terms and inspect only a small annotation/schema artifact. Do not score CARE or promote a sample count until barcode linkage, malignant support, identity provenance, and reference handling are frozen.
