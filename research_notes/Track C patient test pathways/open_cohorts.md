# Track C Published De-identified Patient-Data Audit

Audit date: 2026-10-31. Scope limited to already-published, de-identified GBM/TMZ patient datasets and public cohort manifests. No patient-level records were opened or downloaded.

## Which publicly downloadable patient files could support a real-patient Track C test?

### Takeaway
No publicly downloadable dataset I verified documents the complete record Track C needs: actual TMZ dose and administration dates by course, tumor measurements at treatment decision times, progression or survival, and toxicity. Several open cohorts support narrower imaging or prognosis endpoints under observed care.

### Cited Findings

| Public file(s) | What published source documents | Exact Track C gap |
|---|---|---|
| **UCSD-PTGBM — strongest open endpoint cohort** | Published, de-identified, CC BY 4.0 dataset: 178 GBM patients and 243 post-treatment MRI timepoints; 192 timepoints have residual/recurrent tumor and 51 have treatment changes only. OS and PFS are available for 94 patients. MRI, segmentations, and clinical CSV are distributed through TCIA. [Dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12886825/); [TCIA collection](https://www.cancerimagingarchive.net/collection/ucsd-ptgbm/) | Does not report course-level TMZ administrations/doses, treatment interruptions, or toxicity. Published endpoint labels distinguish tumor from treatment-related change; they do not identify the effect of an adaptive schedule. |
| **ACRIN-DSC-MR-Brain / ACRIN 6677 clinical files** | TCIA lists de-identified clinical XLS/CSV/PDF files as downloadable under CC BY 4.0; 123 trial participants, with 107 having baseline plus at least one post-treatment MRI. Clinical forms include demographics, diagnosis, follow-up, measurement, and time-offset fields. The MRI files are under NIH controlled access. [TCIA data record](https://www.cancerimagingarchive.net/collection/acrin-dsc-mr-brain/) | Published data description does not enumerate patient-level delivered TMZ dates/doses or adverse-event fields in the downloadable subset. Trial publications report fixed bevacizumab+TMZ or bevacizumab+irinotecan assignment and aggregate toxicities. Exact cycle-level drug records remain unverified. |
| **RHUH-GBM** | TCIA lists de-identified clinical CSV and brain-extracted images/segmentations under CC BY 4.0 for 40 GBM patients; original DICOM images require controlled access. It includes preoperative, early postoperative, and recurrence MRI. [TCIA collection](https://www.cancerimagingarchive.net/collection/rhuh-gbm/) | Published manifest does not document actual TMZ course dates/doses or toxicity. Small cohort; no adaptive treatment assignments. |
| **CFB-GBM v4** | Public CC BY 4.0 collection of 264 GBM cases with pre-treatment MRI and post-treatment scans when available. The manifest explicitly lists radiation dose, number of fractions, and delay from baseline MRI to treatment start. [TCIA collection](https://www.cancerimagingarchive.net/collection/cfb-gbm/) | Those treatment fields describe radiotherapy, not actual TMZ dosing by cycle. No toxicity history is listed. |
| **MU-Glioma-Post** | Public CC BY 4.0 clinical XLSX and postoperative MRI/segmentations for 203 glioma patients; source includes treatment and follow-up categories. [TCIA collection](https://www.cancerimagingarchive.net/collection/mu-glioma-post/); [dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12635165/) | Published manifest/paper does not list actual TMZ administration dates/doses, interruptions, or cycle-level toxicity; cohort includes glioma types beyond GBM. |
| **LUMIERE** | Published public dataset of 91 GBM patients, 638 study dates, 2,487 MRI images, expert RANO evaluations, and OS for a subset. Dataset materials state non-commercial use. [Dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9755255/); [file guide](https://springernature.figshare.com/articles/dataset/LUMIERE_dataset_-_Readme_file/21266241) | Treatment course records in the publications do not give actual per-cycle doses, dates, holds, or toxicity. Published standard-care exposure does not test an adaptive policy. |
| **UCSF-ALPTDG** | Public annotated dataset of 298 diffuse-glioma patients, including 80 GBM, with two consecutive follow-up MRIs, expert labels, clinical history, progression consensus, and death interval. Use is subject to its non-commercial data-use agreement. [Dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11294954/) | Paper does not report actual TMZ dose/date by cycle or toxicity. It can support observed-care imaging/prognosis only. |
| **CGGA** | Its public download portal provides clinical/genomic data with therapy and survival information; its image-genomic MRI set is a separate 268-sample resource. [CGGA download page](https://www.cgga.org.cn/download.jsp); [resource paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8498921/) | Published resource descriptions do not specify actual TMZ administration dates/doses, toxicity, or serial decision-time MRI. Therapy/survival fields alone do not reconstruct adaptive dosing. |
| **TCGA/GDC and GLASS public components** | GDC makes most clinical elements openly accessible, with project-specific supplements; GLASS publicly released a resource of standardized somatic variants and clinical annotations, while raw sequencing access is complicated by jurisdiction-specific limits. [GDC clinical-data documentation](https://docs.gdc.cancer.gov/Encyclopedia/pages/Clinical_Data/); [GDC access types](https://gdc.cancer.gov/access-data/data-access-processes-and-tools); [GLASS data-resource page](https://glass-consortium.org/datasets/); [GLASS cohort paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6897368/) | Neither published description supplies actual TMZ delivery by course plus serial tumor measurement and toxicity. GLASS samples are separated surgical specimens, not frequent on/off treatment decision observations. |

### Inferences
- UCSD-PTGBM is the strongest fully open, GBM-specific dataset found for a narrow patient-level MRI endpoint: predict residual/recurrent tumor versus treatment-related change using held-out patients.
- ACRIN 6677 is the closest public clinical-table near miss for a trial-based outcome analysis, but its imaging is controlled and the public manifest does not establish course-level TMZ delivery.
- These data can test imaging classification or prognosis under observed care. They cannot estimate what would happen if a patient followed Track C's alternate dosing rule.

### Gaps
- ACRIN's public data dictionary may contain treatment fields not described in its published manifest. I did not open/download patient files or the patient-level spreadsheets, so actual dose, cycle, and toxicity field availability remains unverified.
- Public dataset manifests do not prove that no related unpublished or request-only records exist.

## Which published de-identified datasets are controlled-access, and could they fill the missing fields?

### Takeaway
The closest controlled-access repeated-imaging resource is SAILOR; the closest controlled-access TMZ trial source is RTOG 0625/ACRIN 6677, with RTOG 0525 as a standard-schedule comparator. None is verified to provide all needed actual exposure fields, and both clinical trials use fixed, protocol-defined strategies rather than tumor-state-guided adaptation.

### Cited Findings
- **SAILOR (controlled, DUA):** 27 high-grade glioma patients with 3–19 serial MRI timepoints, per-timepoint status labels (CRT/TMZ/no/unknown), RANO classes, survival, and days between scans. EBRAINS marks access as controlled and requires a data-use agreement. It does not list TMZ dose-by-cycle or toxicity fields. [EBRAINS record](https://search.kg.ebrains.eu/instances/cae85bcb-8526-442d-b0d8-a866425efff8?noSilentSSO=true)
- **Duke longitudinal MRI (request-only):** 61 GBM patients, 503 MRI sessions, masks, and relative diagnosis/chemoradiation/radiation/surgery/death dates. Files are restricted and request-based; the manifest does not list actual TMZ dose/course dates or toxicity. [Zenodo record and access terms](https://zenodo.org/records/20548528)
- **ACRIN 6677/RTOG 0625 (mixed access):** de-identified clinical tables are public, but the 86.9 GB MRI collection is under NIH controlled access. Trial reports describe recurrent GBM patients randomized to bevacizumab plus dose-dense TMZ or bevacizumab plus irinotecan, MRI at baseline and then about every eight weeks, progression and survival outcomes, and toxicities. The protocol permits toxicity-based TMZ escalation, reduction, delays, or cessation; those are protocol rules, not confirmation that released patient files record each actual dose/date. [TCIA data record](https://www.cancerimagingarchive.net/collection/acrin-dsc-mr-brain/); [trial report](https://pmc.ncbi.nlm.nih.gov/articles/PMC5263144/); [imaging trial paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4383029/)
- **RTOG 0525 / NCTN-NCORP:** trial publication reports 833 patients randomized between fixed standard and dose-dense schedules. NCI archive provides de-identified patient-level data and dictionaries for listed trials through dbGaP; current public listing did not confirm RTOG 0525. NRG requests for eligible data outside the archive require an investigator CV, trial-specific proposal/methods, signed agreement, institutional IRB approval, human-subjects training, and approval. A published secondary analysis requested the number of adjuvant TMZ cycles from RTOG 0525, but does not show that per-cycle administration dates or doses were released. [Trial report](https://pmc.ncbi.nlm.nih.gov/articles/PMC3816958/); [NCI archive](https://dctd.cancer.gov/research/networks/nctn/data-archive); [NRG policy](https://www.nrgoncology.org/wp-content/uploads/2025/09/NRG-Oncology-Data-Sharing-Policy-v3-01Oct19.pdf); [published secondary request](https://pmc.ncbi.nlm.nih.gov/articles/PMC8494241/)
- **EORTC and GBM-X:** EORTC shares published trial data only after a qualified researcher submits a proposal and receives approval. GBM-X contains de-identified IPD from more than 1,200 newly diagnosed patients treated with RT/TMZ, but offers approved investigators analyses rather than direct raw-patient downloads. Neither source's public documentation confirms cycle-level TMZ delivery and toxicity fields. [EORTC access policy](https://www.eortc.org/data-sharing/); [GBM-X paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10681026/)

### Inferences
- SAILOR could support an approved, small-sample analysis of observed MRI change and broad treatment status. Its labels do not identify actual delivered dose or the causal effect of stopping/restarting TMZ.
- ACRIN 6677/RTOG 0625 could support patient-level outcome analysis if the clinical fields are sufficient and the imaging-access request is approved. It cannot validate Track C's policy because treatment was randomized to fixed arms; treatment adjustment was protocol toxicity management, not tumor-state-guided adaptation.
- RTOG 0525 may permit comparison of fixed schedule arms under a qualified investigator. It still cannot test the adaptive rule without the required actions, dose timing, and decision-time covariates.

### Gaps
- No trial-specific data dictionary or case-report form was verified to contain actual delivered TMZ dates/doses, toxicity tied to each course, and the matching imaging/progression timeline.
- A request is not access approval. This project has no confirmed investigator sponsor, IRB determination, or data-use approval for these controlled datasets.
- The current source audit did not establish whether RTOG 0525 is in the NCTN archive or whether NRG would release the required variables.

## What real-patient endpoint can Track C test with open data now?

### Takeaway
The strongest defensible open patient-data study is an imaging endpoint on UCSD-PTGBM: predict residual/recurrent tumor versus treatment-related change on post-treatment MRI, with patient-level holdout. This is real-patient prediction under observed care; it is not an adaptive-TMZ policy test or a claim of treatment benefit.

### Cited Findings
- UCSD-PTGBM is de-identified, publicly available through TCIA under CC BY 4.0, and contains 178 patients/243 MRI timepoints, expert tumor labels, and OS/PFS for 94 patients. [Dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12886825/)
- The labels distinguish 192 timepoints with residual/recurrent tumor from 51 timepoints with post-treatment changes only; OS/PFS coverage is a subset, not all patients. [Dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12886825/)
- ACRIN 6677 clinical data include clinical/outcome forms and progression timing, but its image files require controlled access and the public description does not verify actual course-level TMZ dose records. [TCIA data record](https://www.cancerimagingarchive.net/collection/acrin-dsc-mr-brain/)

### Inferences
- Use UCSD-PTGBM only for a prespecified patient-held-out MRI task, such as classifying tumor versus treatment-related change. Report exactly that endpoint.
- Do not call this TMZ response, treatment-effect estimation, adaptive-policy validation, or patient-specific benefit. Published sources do not establish patient-level actual dosing histories or a policy counterfactual.
- A real adaptive-policy evaluation still needs controlled-access course-level records with actual treatment decisions/doses, matched decision-time imaging, outcomes, and toxicity, or a clinician-led prospective study.

### Gaps
- Patient-level timing of prior TMZ and toxicity are not documented in UCSD-PTGBM's published data record.
- Whether any controlled trial database contains all these fields remains unknown until a qualified investigator reviews its data dictionary and access is approved.

