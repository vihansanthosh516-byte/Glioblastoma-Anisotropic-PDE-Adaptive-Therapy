# Greenwald GBM spatial RNA–protein eligibility gate

**Audit date:** 2026-11-16
**Decision:** **No-go for independent confirmatory validation at current gates; conditional for a small, donor-level exploratory cross-assay check.** No expression matrices, cell-level marker values, or outcomes were read. No large files were downloaded. Existing core/edge benchmark remains a separate endpoint and is not rescued by this candidate.

## What the metadata establishes

- Official Zenodo concept DOI `10.5281/zenodo.8105466` resolves to current version record **12624860**, DOI `10.5281/zenodo.12624860` (updated 2024-08-28). It is an open, published record under **CC BY 4.0**. Attribution is required; reuse must retain the license/citation. GEO GSE237183 is the linked spatial transcriptomics accession.
- Current Zenodo record describes 13 fresh-frozen GBM Visium sections, six IDH-mutant Visium samples, 12 GBM CODEX tissue samples, and one GBM organoid model. The paper reports 13 GBM Visium sections and 12 CODEX sections; it explicitly describes the CODEX sections as **near-adjacent** to Visium. It does not claim same-section or same-cell RNA/protein measurements.
- **Five independent GBM donor identifiers have a donor-level RNA/CODEX overlap** in the published sample/file labels: MGH258, ZH881, ZH916, ZH1007, ZH1019. Other CODEX-labeled donors visible in image names do not have a correspondingly named GBM Visium file in the Zenodo Visium inventory. This is a verified *donor-label overlap*, not an author-verified specimen/section crosswalk. Thus the number currently verified as suitable for donor-level cross-assay validation is **zero**; five are candidates pending the specimen/section map, and five is below the project's six-donor threshold in any case.
- Paper methods describe computational image alignment of CODEX with Visium using STalign. That supports spatial registration as a published method, but does not substitute for proving which CODEX field/section pairs to which Visium section for each donor. Visium coordinates are delivered in per-sample `*_tissue_positions_list.csv` files inside processed archives; CODEX cell segmentation and pseudospot grid are described in the Zenodo README/QuPath package. Coordinates and segmentation appear available, but no inspected public metadata crosswalk binds those coordinate sets to exact matched sections.
- CODEX uses a **40-antibody panel**. The deposit includes `CODEX_antibody_validations.pdf`; this metadata-only audit did not extract the marker list. The README says the QuPath project includes images, IvyGAP pathology annotations, segmentation masks, and a pseudospot grid. RNA is spot-based Visium (55 μm capture spots; mixed cells per spot), so this is tissue/region-scale orthogonal validation, not single-cell RNA–protein concordance.

## Access and burden

Access class is **A — direct public download** for the Zenodo files; no account or case-by-case approval is indicated. Zenodo API record `https://zenodo.org/api/records/12624860` lists **38 files totaling 25,333,309,514 bytes (~25.3 GB)**. The large image/project payloads include `CODEX_qupath_project_v2.zip` (**2,219,707,544 bytes**), `MGH258_v6.ome.tif` (**4,468,568,903 bytes**), and multiple other OME-TIFFs around 0.7–2.6 GB. No bulk download is recommended at this gate.

Small metadata files inspected directly (10,730-byte `visium_dataset_description.csv`; 556-byte `readme_CODEX.txt`):

- `https://zenodo.org/api/records/12624860/files/visium_dataset_description.csv/content` — identifies Visium bundle names and processed contents; bundles include `GBM_MGH258.tar.gz` (25,901,277 bytes), `GBM_ZH881T1.tar.gz` (46,130,044), `GBM_ZH916T1.tar.gz` (49,897,943), `GBM_ZH1007nec.tar.gz` (22,248,502), and `GBM_ZH1019T1.tar.gz` (23,447,164), among other regions/sections.
- `https://zenodo.org/api/records/12624860/files/readme_CODEX.txt/content` — says the compressed CODEX archive contains processed CODEX images, corresponding H&E scans, and a QuPath project with images, pathology annotations, segmentation masks, and pseudospot grid.
- Full API listing with exact paths, bytes, MD5 values, and per-file URLs: [`spatial_manifest.json`](../output/paper_strengthening_metadata/spatial_manifest.json).

The Zenodo record is the controlling stated license for deposited files. GSE237183’s GEO record links the RNA accession but does not itself establish the RNA–CODEX specimen map. No intended-use conflict was found in the inspected deposit metadata; this is not a legal determination.

## Gate and unresolved mapping

The paper supports **near-adjacent-section** RNA/protein comparison, not same-section matching. Shared donor codes make a donor-level join plausible for five tumors, but the inspected public README, Visium description CSV, GEO record, and paper do not provide a complete table mapping, for each donor, the exact Visium sample/region to the CODEX file, slide, and image/field. Therefore:

1. **Same section:** not established; source paper says near-adjacent sections.
2. **Exact adjacent-section/sample pairing by donor and region:** not verified in public metadata inspected.
3. **Eligible independent donors:** 0 verified; 5 donor-code-overlap candidates, all dependent on the unresolved section map. Even if all five map, this is below six donors.
4. **Independent validation:** not independent of the paper that generated and analyzed these data; the proposed use would be a reanalysis/orthogonal endpoint on the same study cohort.

Proceed only if the public supplement/code supplies an unambiguous per-donor section/file crosswalk and a new protocol specifies its own endpoint, precision plan, and donor requirement. The six-donor rule belongs to the existing regional benchmark; it is not a universal threshold for every new biological question. Five candidate donors nevertheless give limited independent replication, and multiple sections cannot increase that donor count. Otherwise stop. Do not describe this as confirming the separate core/edge benchmark.

## Sources

- Greenwald et al., *Cell* (2024), [PMC11088502](https://pmc.ncbi.nlm.nih.gov/articles/PMC11088502/), especially Fig. 1 and data-availability statement.
- [GEO GSE237183](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237183).
- [Zenodo concept record, DOI 10.5281/zenodo.8105466](https://zenodo.org/records/8105466); [current version record 12624860](https://zenodo.org/records/12624860); [official API metadata/file listing](https://zenodo.org/api/records/12624860).
