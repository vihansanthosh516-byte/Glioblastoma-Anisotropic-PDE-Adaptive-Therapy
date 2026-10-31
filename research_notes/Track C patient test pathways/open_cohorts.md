# Track C Real-Patient Data Pathways

## Which named longitudinal glioblastoma cohorts contain the fields needed to test the adaptive TMZ policy?

### Takeaway
No public dataset manifest or published cohort description I checked reports the full combination needed for a patient-level test of Track C: actual time-stamped TMZ administrations and dose changes, repeated decision-time tumor measurements, and outcome plus toxicity. Several cohorts can support real-patient imaging or prognosis analyses under observed care, but those analyses cannot establish benefit from the adaptive rule.

### Cited Findings

| Cohort | Patient-level evidence reported | Access and exact gap for Track C |
|---|---|---|
| **UCSF-ALPTDG** | 298 diffuse-glioma patients, including 80 GBM; two consecutive follow-up MRIs per patient, expert voxelwise labels, clinical history, progression consensus, and death interval. | Public noncommercial data use is subject to a data use agreement. The paper does not describe actual per-cycle TMZ dates, delivered doses, holds, or dose reductions. Suitable for observed-care imaging-change/prognosis work, not adaptive-policy effects. [Dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11294954/) |
| **LUMIERE** | Longitudinal GBM MRI with RANO evaluations and clinical information; published cohort describes standard Stupp treatment. | Published materials reviewed do not specify a per-patient administration calendar, actual dose received by cycle, treatment interruptions, or toxicity linked to each decision time. Reuse terms need direct resolution before analysis. [Dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9755255/) |
| **SAILOR** | EBRAINS record describes 27 high-grade-glioma patients, 3–19 MRI time points, survival, RANO classes, and per-scan treatment-status labels such as chemoradiation, TMZ, no treatment, or unknown. | This gives coarse treatment status, not dated actual dose/course records or toxicity. The accessible metadata did not establish exact file-level access conditions or reuse license. [EBRAINS record](https://search.kg.ebrains.eu/instances/cae85bcb-8526-442d-b0d8-a866425efff8?noSilentSSO=true) |
| **TCGA/GDC** | GDC supports patient clinical exposures and project-specific clinical supplements; TCGA supplements can be in different formats. | The common clinical record does not establish complete actual TMZ delivery by cycle, treatment holds, or repeated decision-time imaging. A cohort-specific field audit would be needed; TCGA/GDC alone is not an adaptive-policy dataset. [GDC clinical-data documentation](https://docs.gdc.cancer.gov/Encyclopedia/pages/Clinical_Data/) |
| **GLASS** | Longitudinal molecular samples and clinical annotations from tumors collected at separate surgeries; public GLASS resource reports multi-timepoint profiles from more than 400 glioma patients. | Repeated tumor samples are sparse surgical time points, not cycle-level TMZ exposure or serial decision-time burden. The public resource does not establish actual drug administration, toxicity, or adaptive treatment assignments. [GLASS resource](https://glass-consortium.org/); [longitudinal cohort paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6897368/) |
| **MU-Glioma-Post** | Public TCIA collection: 203 subjects, postoperative MRI/segmentations, clinical XLSX; collection lists CC BY 4.0. | Clinical file supports diagnosis, treatment, and follow-up fields, but source descriptions do not document actual per-cycle TMZ administration dates/doses, holds, or dose reductions. It can support retrospective imaging/progression work under observed treatment, not a policy-effect estimate. [TCIA collection](https://www.cancerimagingarchive.net/collection/mu-glioma-post/); [dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12635165/) |
| **Duke longitudinal GBM MRI (2026)** | 61 GBM patients, 503 MRI sessions, 3–18 sessions each, tumor masks, scan dates, and relative dates for diagnosis, chemoradiation, radiation, surgery, and death. | Files are restricted and require a request. Description does not list per-cycle TMZ administration, actual dose, interruptions, or toxicity. It is a strong patient imaging/prognosis resource, not an adaptive TMZ test. [Zenodo record and terms](https://zenodo.org/records/20548528) |
| **CFB-GBM v4** | TCIA collection reports 264 GBM cases with pre/post-treatment imaging and radiotherapy details, including radiation dose/fractions and treatment-initiation delay. | Described treatment fields concern radiotherapy; the manifest does not report actual TMZ dose/course timing. Public imaging can test imaging endpoints under observed care, not TMZ policy benefit. [TCIA collection](https://www.cancerimagingarchive.net/collection/cfb-gbm/) |

### Inferences
- A real-patient study using these imaging cohorts could test whether a frozen model predicts later MRI change, progression, or survival under the treatment patients actually received. That is prognostic validation, not a test of what would have happened under Track C's alternative dosing decisions.
- Without time-stamped treatment choices and delivered doses, the dataset cannot reconstruct whether a patient was on treatment or off treatment at each policy decision, or estimate effects of following the adaptive rule.
- Even complete observational dosing histories would need a prespecified causal design and enough treatment variation; a prospective policy trial would be needed to test actual clinical use and safety.

### Gaps
- Published descriptions and manifests do not prove that no unpublished or request-only dose records exist. A data dictionary and case-report form must be obtained before selecting a cohort.
- I did not access or download patient-level data. File-level dates, dose fields, toxicity fields, missingness, de-identification, and permitted school use remain unverified.
- The Duke files are request-only, and the project records that a request was submitted on 2026-10-28; access and school-use approval remain unconfirmed.

## Is there a stronger requestable patient cohort or randomized comparator?

### Takeaway
The strongest clinical-data route found is an investigator-mediated request for randomized TMZ-trial data, especially RTOG 0525 if its patient-level records contain actual delivered courses. It compares fixed schedules, so it cannot validate Track C's adaptive policy by itself; current access and schedule-detail availability are unconfirmed.

### Cited Findings
- RTOG 0525 randomized 833 newly diagnosed GBM patients to fixed standard or dose-dense TMZ schedules. It reported overall survival and progression-free survival, and higher grade 3+ toxicity with dose-dense treatment; it did not assign a tumor-burden-guided on/off policy. [Primary trial report](https://pmc.ncbi.nlm.nih.gov/articles/PMC3816958/)
- NCI's NCTN/NCORP Data Archive provides de-identified patient-level clinical data, data dictionaries, and variables used in published analyses; some studies also have imaging or genomic data in dbGaP. The current list includes two primary CNS-tumor trials, but the reviewed material did not confirm that RTOG 0525 is among them. [NCI archive](https://dctd.cancer.gov/research/networks/nctn/data-archive)
- NRG's policy directs archive-deposited trials to NCI. For eligible NRG trials outside the archive, a request requires an investigator CV/biosketch, trial-specific aims and analysis methods, signed agreement, institutional IRB approval, human-subjects training, and panel approval. [NRG data-sharing policy](https://www.nrgoncology.org/wp-content/uploads/2025/09/NRG-Oncology-Data-Sharing-Policy-v3-01Oct19.pdf)
- Ipi-Glio's primary report states that de-identified individual data and a data dictionary may be requested through Oxford OCTO with a reason, investigator support, proposal approval, and a data-access agreement. The trial used standard fixed adjuvant TMZ schedules and MRI follow-up; the publication does not establish that released records include actual administered dose/date histories needed for an adaptive-policy emulation. [Primary trial report](https://pmc.ncbi.nlm.nih.gov/articles/PMC12342986/)
- GBM-X pools individual-level data from more than 1,200 patients in six trials and two institutional cohorts, but its paper describes investigator-requested analyses rather than an openly downloadable patient-level dataset. Its cohorts received standard radiotherapy/TMZ. [GBM-X paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10681026/)

### Inferences
- RTOG 0525 could support a real-patient comparison of fixed schedule arms and might provide a cohort for prognosis under randomized treatment. It cannot directly answer whether an adaptive policy improves outcomes because the policy was never assigned.
- A qualified neuro-oncology investigator could first request the RTOG 0525 data dictionary and confirm actual TMZ start/stop dates, each delivered dose, treatment holds/reductions, imaging dates/progression, survival, and toxicity. Only then can the team decide whether a valid patient-data analysis is possible.
- A real intervention test of an adaptive schedule requires clinician-led prospective research, IRB oversight, and patient consent. A retrospective dataset cannot establish the safety of changing a patient's TMZ schedule.

### Gaps
- RTOG 0525's presence in the currently listed NCTN/NCORP archive was not verified; if absent, availability through NRG is not guaranteed.
- No public data dictionary or case-report form located here confirms actual delivered course-level TMZ exposure, rather than assigned regimen or intended schedule.
- No direct access request was submitted. The project's current status does not include an approved institutional investigator, IRB determination, or cohort data-use approval.

## What patient-level claim can Track C support now, and what exact gate comes next?

### Takeaway
Track C cannot claim patient-specific benefit or validate adaptive dosing from the cohorts described. The concrete next gate is a clinician/investigator obtaining and auditing a randomized GBM TMZ-trial dictionary for delivered-dose and outcome fields; if those fields are absent, use patient data only for a clearly labeled prognosis study or pursue a prospective trial.

### Cited Findings
- The project's current plan labels Track C outputs as simulation using assumed kill rates, not a real-drug evaluation, and records that model testing stops pending measured transition data and review. [Project START_HERE](../../START_HERE.md); [Track C sensitivity report](../../paper/trackc_sensitivity_screen_2026-10-30.md)
- The NCI archive describes its files as de-identified trial clinical data and data dictionaries, while NRG's outside-archive process requires investigator and institutional review documentation. [NCI archive](https://dctd.cancer.gov/research/networks/nctn/data-archive); [NRG data-sharing policy](https://www.nrgoncology.org/wp-content/uploads/2025/09/NRG-Oncology-Data-Sharing-Policy-v3-01Oct19.pdf)

### Inferences
- The best immediate patient-data test is not an adaptive-treatment claim. It is a frozen, patient-held-out test of an imaging/progression predictor on an eligible cohort, clearly named as prediction under observed treatment.
- To test policy value, the necessary record must include each treatment decision and actual delivered exposure aligned to scan dates and progression/toxicity, with treatment variation sufficient for a defensible causal comparison; otherwise only a clinician-led prospective study can test the policy.
- Do not open another sealed simulation test to compensate for missing clinical data. Patient data access and field sufficiency are separate gates from model performance.

### Gaps
- No cohort reviewed here is confirmed to contain the complete exposure/outcome record required for Track C's policy estimand.
- Any cohort access must be checked for de-identification, license/data-use terms, institutional/school approval, and whether those terms allow this specific analysis before patient-level use.
