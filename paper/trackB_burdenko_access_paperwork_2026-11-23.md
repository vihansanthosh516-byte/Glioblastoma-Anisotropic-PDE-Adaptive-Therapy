# Burdenko GBM Progression: access paperwork (draft, 2026-11-23)

Status: DRAFT. Placeholders in [BRACKETS] must be filled by the applicant. Nothing here is sent.

## 1. Rules found (source: https://www.cancerimagingarchive.net/nih-controlled-data-access-policy/, read 2026-11-23)
- TCIA no longer hosts controlled data. Access is requested in dbGaP. Study: "Radiology Imaging Data Containing Faces from TCIA", **phs004225**.
- The requester must be tenure-track professor level or senior scientist. Students, trainees and lab staff may NOT submit the request. They may work on the project if an eligible person oversees it.
- Researchers at institutions in China, Russia, Iran, North Korea, Cuba and Venezuela are barred (NOT-OD-25-083).
- After approval: NCI Data Commons login, make API key, get manifest from TCIA Dataset Explorer, download with TCIA Data Retriever.
- Dataset page: https://www.cancerimagingarchive.net/collection/burdenko-gbm-progression/ (180 patients).

## 2. Unverified (check before sending)
- dbGaP usually also needs an institutional Signing Official. Not confirmed on the TCIA page.
- Review time is unknown.
- Whether the images give 2 or more timed follow-up scans per patient is not matched to our clinical table (657 dates, median 3).
- Whether a school (not a hospital or university) can be the applicant institution.

## 3. Email to the applicant (adult)
Subject: Request: faculty applicant for NIH brain MRI data access (research project)

Hello [NAME],
I am doing a computational research project. It asks: from a post-treatment brain MRI plus scan dates and treatment dates, how big will the tumor core be at the next scan? I need the Burdenko glioblastoma MRI images. NIH rules say the access request must come from a faculty-level person. Could you be the applicant? I will write the text, run the analysis and keep the data safe. I need about [30 minutes] of your time to review and sign. The dataset is here: https://www.cancerimagingarchive.net/collection/burdenko-gbm-progression/
Thank you, [YOUR NAME]

## 4. Research Use Statement (draft, about 250 words)
Title: Predicting next-scan tumor-core volume after glioblastoma treatment from serial MRI.

Objective: Test whether a method that uses a post-treatment MRI, the scan dates and the treatment dates predicts the log tumor-core volume at the next scan better than simple baselines (carry the last volume forward, population growth rate, exponential and Gompertz growth fits, gradient boosting).

Data use: We will use only the Burdenko images (T1 with contrast, FLAIR) and the open clinical table with scan dates. A segmentation model that is frozen and already trained on other public data will outline the tumor. No model is trained on these images. The analysis plan, the method and the scoring rule are written and dated before the images are opened. The images are scored once.

Why this data: it has dated follow-up scans after radiotherapy. Most public sets have no dates.

Non-identification: We will not try to identify any person or reconstruct faces. We will not link the data to other data that could identify a person. Images stay on [SECURE COMPUTER / INSTITUTION SERVER]. We will not share the images. We will share code and summary results only.

Outputs: a report and code. No clinical use, no treatment advice.

## 5. Data security plan (draft)
- Storage: encrypted disk, [LOCATION], password protected, access by named people only: [APPLICANT], [YOUR NAME].
- No cloud copy of raw images unless the institution allows it and the data rules permit it. (We used a private Kaggle dataset for earlier cohorts. Do NOT upload Burdenko images there unless dbGaP terms allow it. Check first.)
- Delete the images at project end or when the access term ends.

## 6. What the student does vs what the applicant does
- Student: writes text, runs analysis, keeps the data safe.
- Applicant: owns and submits the dbGaP request, signs, is responsible under the terms. Only they can.

## 7. Next step after approval
Freeze VIG unchanged. Pre-register. Open the cohort once. Report with the 4-question format.
