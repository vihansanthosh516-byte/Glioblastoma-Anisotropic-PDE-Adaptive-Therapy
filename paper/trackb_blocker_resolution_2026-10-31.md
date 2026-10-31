# Track B blocker audit and resolution

**Date:** 2026-10-31  
**Scope:** Existing Track B PDE tumor-growth forecasting core. CFB annotation-change remains a separate exploratory side study.  
**Access rule:** No author or dataset-owner requests; do not use controlled files.

This audit combined a repository/protocol review with two independent GPT-6 Luna source screens: one for CFB evidence gates and one for direct-public PDE cohorts. No large data were downloaded, no access request was submitted, and no patient score was run.

## Decision

No new dataset download is needed. MU-Glioma-Post is already the best-fitting direct-public cohort for the existing PDE workflow, and its data are already in the project. CFB-GBM v4 is also downloaded locally, but its selected follow-up GTV masks are generated labels and its unresolved evidence gates remain open. Public alternatives do not provide a cleaner, untouched, direct-public replacement for patient-level PDE validation.

The core MU run is **not cleared for confirmatory patient scoring**. The available code and synthetic validations are ready; the remaining blockers concern source evidence and patient-level validation, not a missing download or a missing model implementation. Do not turn a gate into a pass by changing its cutoff after seeing outcomes. Existing MU and LUMIERE results remain exploratory because those cohorts have already been examined.

## Gate ledger

| Gate | Current finding | Status |
|---|---|---|
| Dataset availability and rights | MU-Glioma-Post is publicly downloadable through TCIA; the project already has the cohort. TCIA lists the image/segmentation package and clinical files under CC BY 4.0. A free TCIA account is required. | **Closed for access** |
| Candidate feasibility | The local metadata screen found 80 exact-GBM patients with at least three mask-available visits; 18 had the planned 76–104-day second-to-third interval. Headers and MRI-only screens exist. | **Closed for feasibility only** |
| Untouched validation | MU, LUMIERE and UCSD have already appeared in earlier Track B analyses. They cannot be presented as new independent confirmation. | **Open; no eligible replacement found** |
| Registration | MU source material describes rigid atlas registration and general image-quality review, but not case-level longitudinal QC, selected-patient transform records, or tumor-boundary alignment. CFB likewise has a cohort-level registration method but no per-case acceptance record. Matching headers do not prove anatomy is aligned. | **Open** |
| Target support and negative meaning | MU's source defines labels 1–4 but does not define label 0 as a valid negative over the whole image grid. The existing MRI nonzero mask is explicitly only a support proxy. CFB's brain mask and matching headers do not establish annotation coverage or valid negatives. | **Open** |
| Label provenance | MU masks are automated/derived with reported refinement; that does not make them independent truth. In the selected CFB candidates, all t1/t2 GTVs are `genSegmentation`; the paper reports 44 reviewed generated masks but gives no file IDs, so reviewer membership is unknown. | **Open for biological validation** |
| Timing and treatment | MU offers day-valued imaging intervals. CFB's selected 38-case interval values come from week-resolution metadata (77/84/91/98 days), not verified exact scan days. Treatment changes can affect observed growth and must remain part of the interpretation. | **MU interval data closed; CFB exact-day gate open** |
| Scoring implementation | The gated MU loader, patient-level runner, comparators, fixed-horizon option, support-aware scoring, and synthetic checks already exist. The latest metadata-only preflight stopped before opening target masks because registration/support evidence was missing. | **Implementation ready; patient scoring blocked** |

## Dataset screen

- **MU-Glioma-Post — keep for Track B development.** It is the closest fit to the existing postoperative PDE question and is already local. The project has prior patient-level analyses, so future work on it is exploratory/development evidence, not an untouched test.
- **CFB-GBM — keep local; do not redownload.** It has three planned visits and treatment metadata, but 76/76 selected future labels are generated, the per-mask reviewer subset is unidentified, exact t1-to-t2 days are unavailable in the joined table, and registration/support evidence is incomplete. Its binary GTV is enhancing plus necrotic tumor and excludes edema; that must not be silently equated with a broader PDE target.
- **RHUH-GBM — do not download for this PDE aim.** Its public NIfTI release has expert-corrected labels, but the three visits are preoperative, within 72 hours after surgery, and at recurrence. The first transition is dominated by resection, which is incompatible with smooth untreated PDE growth fitting.
- **LUMIERE — do not redownload.** Its CC0 release has many irregular follow-up scans and automated masks, but this cohort has already been used in the project and cannot provide fresh confirmation.
- **UCSD-PTGBM — do not treat as fresh confirmation.** It is public, but it was already used for the separate UCSD annotation test and is not an untouched Track B cohort.
- **Duke, SAILOR, and Burdenko — exclude under the no-request rule.** The relevant files require restricted access, a data-use agreement, or a signed license request. No request will be sent and no restricted files will be downloaded.

## What this audit resolves

1. **Download uncertainty is closed:** keep using the local MU and CFB data; do not download Duke, RHUH, or another dataset for this aim.
2. **The no-contact constraint is applied:** the older MU draft email is not the next step. Public-source checks found no existing case-level evidence that closes the MU registration or support gates.
3. **The CFB review preparation is complete, but the review is not:** all 38 candidate sets have 5-mm navigation pages and a blank pair-level form. No independent full-resolution review has occurred, and the pages do not count as a registration pass.
4. **No patient result is being manufactured:** no CFB GTV/brain-mask voxel values or new patient forecasts were accessed in this audit. Prior patient results retain their existing exploratory status.

## Three workstreams completed without scoring

### 1. MU-Glioma-Post readiness

The data are present in the main checkout at `C:\Users\vihan\20206 science fair\data\tcia\MU-Glioma-Post`, not in this task worktree. A filename-only inventory found 2,978 compressed NIfTI files totaling 11,890,059,719 bytes. Readback of the existing metadata and MRI-only reports confirmed 18 selected 76–104-day candidates, 72 headers checked, 216/216 MRI sequence files present for the selected cohort, 54 T1n MRI volumes previously read for image-only QC, and zero tumor-mask files or voxel arrays opened by those two screens. The MRI-only report explicitly remains descriptive: registration and valid-support gates are not cleared. This worktree initially had Git LFS pointers for the two small JSON reports; only those reports were hydrated, and both now parse as valid JSON. The patient data remain in the main checkout, so a future run from this worktree must point to that data root explicitly.

### 2. CFB-GBM readiness

The v4 package directory at `C:\Users\vihan\Downloads\PKG - CFB-GBM version 4\CFB-GBM` contains 3,976 `.nii.gz` files totaling 208,424,306,322 bytes. The separate `C:\Users\vihan\Downloads\CFB-GBM` directory has only 50 completed NIfTI files and contains `.partial`/Aspera checkpoint files; use the v4 package directory for the existing candidate work and do not restart the download. The pre-existing header audit read 342 headers for 38 candidates: 38/38 shared grids, zero header errors, and no GTV or brain-mask voxel reads. The 760-page MRI aid is complete, but no independent full-resolution reviewer has assessed it. TCIA's reported 3,981 images versus the local package's 3,976 NIfTI files remains an aggregate discrepancy; all selected candidate paths were present in the earlier filename-only check.

### 3. Alternative-dataset disposition

Local inventory confirms RHUH-GBM NIfTI is already present (720 files; 3,098,204,081 bytes), as are UCSD-PTGBM (4,047 files; 35,582,320,168 bytes) and the UCSD BraTS 2024 test package (1,322 files; 11,883,676,648 bytes). RHUH's three timepoints still do not fit smooth PDE fitting across the surgery-dominated preoperative-to-early-postoperative transition. UCSD has already been used for the separate annotation test. LUMIERE has already been examined in Track B. Duke, SAILOR, and Burdenko remain excluded under the no-request rule. Therefore no alternate package should be downloaded for this aim.

## Forward path without contacting data owners

1. Preserve the locked Track B PDE question and the existing MU implementation. Report prior patient results as exploratory and keep the confirmatory claim closed because no untouched direct-public cohort passed the screen.
2. If additional MU patient scoring is still needed for development, first amend the protocol to state that the target is the **released segmentation label**, not biological tumor truth, and define an evaluation support restricted to a source-justified brain domain. Keep origin-only PDE support separate. Do not infer this support from foreground voxels or tune it on outcomes.
3. Freeze a target-blind registration QC procedure for the selected MU cases. Use only MRI/anatomical evidence, report every case and failure, and do not call atlas-space headers or similarity summaries a pass by themselves.
4. Only after those two protocol decisions are documented may the already-built runner be considered for an explicitly exploratory MU run. Any result will be reported against the released masks, alongside persistence and geometric comparators, with the cohort reuse and annotation limits visible.
5. Keep CFB out of PDE scoring unless public source evidence or a predeclared, independent MRI-only review closes its registration and support gates. Do not change the endpoint to obtain a positive result.

## Primary sources checked

- [MU-Glioma-Post TCIA record](https://www.cancerimagingarchive.net/collection/mu-glioma-post/) and [dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12635165/).
- [CFB-GBM TCIA record](https://www.cancerimagingarchive.net/collection/cfb-gbm/) and [CFB-GBM v2.0 paper](https://arxiv.org/html/2608.17884).
- [RHUH-GBM TCIA record](https://www.cancerimagingarchive.net/collection/rhuh-gbm/) and [dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10551826/).
- [LUMIERE Figshare release](https://springernature.figshare.com/articles/dataset/LUMIERE_dataset_-_MRI_data_and_automated_segmentations/21249516) and [dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9755255/).
- [UCSD-PTGBM TCIA record](https://www.cancerimagingarchive.net/collection/ucsd-ptgbm/).
- [Duke longitudinal GBM record](https://zenodo.org/records/20548528), [SAILOR record](https://search.kg.ebrains.eu/instances/cae85bcb-8526-442d-b0d8-a866425efff8?noSilentSSO=true), and [Burdenko-GBM TCIA record](https://www.cancerimagingarchive.net/collection/burdenko-gbm-progression/).
