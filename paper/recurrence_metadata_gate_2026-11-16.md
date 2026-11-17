# GSE174554 recurrence metadata gate — 2026-11-16

## Decision
**Conditional; do not protocol-lock a fixed pair count yet.** Official GEO metadata confirms the series advertises 80 human IDH-wild-type GBM snRNA-seq specimens (40 primary and 40 matched recurrent). The sample-level GEO listing is broader and its pair labels contain duplicates and incomplete pairs. Resolve the sample crosswalk and the 36-versus-40 scope before selecting analysis units. No expression matrices, patient outcomes, or large files were opened or downloaded.

## Verified from primary sources and official GEO metadata

- Wang et al., “A single-cell atlas of glioblastoma evolution under therapy reveals cell-intrinsic and cell-extrinsic therapeutic targets,” *Nature Cancer* (2022), [paper](https://www.nature.com/articles/s43018-022-00475-x), [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9767870/).
- The paper reports 86 specimens from 49 patients overall and describes 36 patients with paired untreated primary and matched first-recurrence specimens. It separately notes four unmatched untreated-primary and six unmatched first-recurrence specimens. These are source-paper cohort/analysis scopes; do not force them to reconcile arithmetically with current GEO listings.
- Official GEO series metadata for [GSE174554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174554) says it profiled **80 human IDH-wild-type GBM specimens by snRNA-seq: 40 primary and 40 patient-matched recurrent**. Series status is public; last update shown is 2023-02-02.
- The retrieved official SOFT series record enumerates **113 GSMs**: **111 human and 2 mouse** by sample organism metadata. **81 sample titles contain `_snRNA`**. The GEO series’ 80-specimen description therefore must not be assumed to mean every `_snRNA` title is one unique eligible specimen.
- Pair-manifest audit of title/characteristics fields identifies **36 numeric pair IDs** among snRNA-titled samples. Of these, **27 IDs have exactly one sample labeled Primary and one labeled Recurrent** in the parsed metadata. This is a literal metadata pattern, not a certified usable-pair count: IDs #10 and #11 each also have a recurrent `v2` library, ID #39 has one primary and two recurrence specimens with different SF identifiers, six additional IDs are recurrence-only, and 12 snRNA-titled samples have `pair#: NA`. Duplicate libraries/specimens and first-recurrence selection need source-level resolution before eligibility is assigned.
- Thus, the original paper’s 36 analyzed patient pairs, GEO’s series-level 40 advertised pairs, and this sample-metadata parse’s 27 unambiguous one-primary/one-recurrence IDs are distinct scopes. **Do not substitute 27 for the final usable count, and do not treat the advertised 40 as independently verified eligible pairs.**
- GEO series metadata lists supplementary paths: `GSE174554_RAW.tar`, `GSE174554_Spatial_proteomics_matrix.txt.gz`, `GSE174554_Spatial_transcriptomic_matrix.txt.gz`, and `GSE174554_Tumor_normal_metadata.txt.gz` under `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174554/suppl/`. File sizes and file-level terms were not present in the metadata inspected and remain unverified. The series record also states raw data are unavailable due to privacy concerns, while the paper’s data-availability statement says raw sequenced reads are in EGA [EGAS00001004909]; actual EGA access terms were not audited. GEO raw-tar naming should not be taken as proof that raw reads are publicly downloadable.
- Assay is single-nucleus RNA sequencing. GEO describes selected cohorts with scATAC-seq, digital spatial profiling, and spatial transcriptomics. The sample-level overlap of these modalities with paired RNA specimens was not fully enumerated here.
- The paper defines GBM as IDH-wild-type GBM and says four IDH1-mutant tumors were excluded from GBM analyses. This is not a per-sample molecular-status crosswalk. GEO series-level IDH-wild-type wording also does not prove each GSM’s IDH status.

## Donor, timepoint, assay, and origin limits

GEO sample metadata provides GSM, title/sample identifiers (including SF IDs), progression labels, diagnosis, age, gender, pair number, platform, organism, and library strategy. The pair number is a submitter-supplied metadata field, not independently verified biological identity. Inconsistencies in labels (e.g., `GBM` versus `Brain tumor`, capitalization, and `v2` titles) mean the donor/timepoint/assay mapping needs reconciliation against Supplementary Table 1 and the original analysis before eligibility is finalized. GSE174554 is the Wang cohort itself: use it as a reanalysis/robustness cohort, not independent replication of Wang’s reported biology. Project metadata records it as a proposed discovery cohort but does not establish patient-level overlap with other project cohorts.

## Gate to close before protocol lock

Compare official GSM/SOFT records with Supplementary Table 1 and the paper’s analyzed set. Resolve donor identity and recurrence specimen selection, distinguish technical `v2` libraries from specimens, reconcile IDH status, and document the exact processed file(s), sizes, and terms. Keep raw-read access separate (EGA) from GEO processed/supplementary files. This is still metadata-only work.

## Local audit artifacts

- Official supplementary-directory listing was also retrieved without opening data files: [recurrence_file_listing.html](../output/paper_strengthening_metadata/recurrence_file_listing.html). Its rounded sizes are 1.5G for `RAW.tar`, 6.6K for spatial proteomics, 472K for spatial transcriptomics, and 1.4M for tumor/normal metadata. These are server display sizes, not exact byte measurements. Archive naming does not establish its contents, and no archive or expression matrix was opened.
- Official series SOFT text: [recurrence_series_metadata.txt](../output/paper_strengthening_metadata/recurrence_series_metadata.txt)
- Sample-level parsed GEO manifest: [recurrence_sample_manifest.json](../output/paper_strengthening_metadata/recurrence_sample_manifest.json)
- Pair-ID grouping and caveats: [recurrence_pair_manifest.json](../output/paper_strengthening_metadata/recurrence_pair_manifest.json)

## Sources and retrieval limits

- [Wang et al. (2022), Nature Cancer](https://www.nature.com/articles/s43018-022-00475-x)
- [Wang et al. (2022), PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9767870/)
- [NCBI GEO GSE174554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174554) — official metadata retrieval is preserved in the local SOFT artifact above; GEO website rendering may CAPTCHA.
- [EGA EGAS00001004909](https://ega-archive.org/studies/EGAS00001004909) — paper-cited location for raw reads; access terms not audited.

No expression matrices or outcomes were read. Counts above describe reported claims and explicit metadata parsing, not inferred eligible patients.
