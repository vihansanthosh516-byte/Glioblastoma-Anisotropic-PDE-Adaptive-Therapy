# Data-request drafts (2026-11-12)

Drafts only. Nothing here was sent. The author sends each message after review. Fill every `[bracket]`. Plan and evidence: `paper/data_request_plan_2026-11-12.md`.

Rules for all three: say you are a student; do not overstate the project; ask for the smallest useful thing; accept any terms they set; save every reply and every agreement in `paper/`.

---

## Draft 1: Moniz Garcia et al. 2026 (Mayo Clinic Florida)

Send to: the corresponding authors listed on the paper (https://pmc.ncbi.nlm.nih.gov/articles/PMC13550674/).

**Subject:** Student request: de-identified data from your GBM multi-sampling paper (Neuro-Oncology 2026)

Dear Dr. [name],

My name is [name]. I am a [grade] student at [school] in [city, country]. I am doing an independent computational research project on glioblastoma progression, [with / without] a mentor ([mentor name and institution, if any]).

I read your paper "Multi-sampling allows intra-tumoral heterogeneity querying and vulnerability profiling in glioblastoma". It is the only study I found that has MRI-guided regional samples, molecular profiles, and measured drug tests for the same patients. My project has these three data types only in separate patient groups, so I cannot link them.

Your data statement says sequencing data are available on request. May I ask if the following could be shared in de-identified form?

1. The gene-level expression table for each regional sample, with its region label (core or margin).
2. The per-sample drug-response values for temozolomide.
3. If it exists, the tumor volume or status at any follow-up MRI. I do not need the images.

I would use the data only to describe, case by case, how regional expression relates to measured temozolomide response. With six patients I would report it as a descriptive case series, not as a statistical test or a clinical claim. I would cite your paper, follow any terms you set, share nothing further, and send you my write-up before I show it to anyone.

I understand that institutional rules may not allow this. If so, I would be grateful for any advice.

Thank you for your time.

[name]
[school]
[contact]

---

## Draft 2: Ntafoulis et al. 2023 (Erasmus MC)

Send to: the corresponding author listed on the paper (https://www.nature.com/articles/s41416-023-02402-y).

**Subject:** Student request: per-culture temozolomide values for GSE232173

Dear Dr. [name],

My name is [name]. I am a [grade] student at [school] in [city, country]. I am doing an independent computational research project on glioblastoma.

I read your paper "Ex vivo drug sensitivity screening predicts response to temozolomide in glioblastoma patients and identifies candidate biomarkers". Thank you for making the RNA-seq counts public as GSE232173.

May I ask two things?

1. Is there a table of the measured temozolomide response (AUC or percent viability at 100 µM) for each of the 24 cultures in GSE232173, matched by culture ID? If it is already in the supplement, please point me to the file.
2. Your data statement mentions further RNA-seq data on reasonable request. Could the cultures with both RNA-seq and a temozolomide measurement be shared in de-identified form?

I would use the data only to test, with a protocol written in advance, whether an expression score is associated with measured in-vitro temozolomide response. I would report it as exploratory and in-vitro only, with no patient or resistance claim. I would cite your paper, follow any terms you set, and share nothing further.

Thank you for your time.

[name]
[school]
[contact]

---

## Draft 3: Request to a local professor for sponsorship

Send to: a faculty member in neuro-oncology, radiology, medical physics, or biomedical data science at [nearby university or hospital]. Pick someone who publishes with public imaging data.

**Subject:** Student seeking a faculty sponsor for a dbGaP data request (glioblastoma MRI)

Dear Professor [name],

My name is [name]. I am a [grade] student at [school]. For about [time] I have worked on an independent computational project on glioblastoma progression using public data.

I have reached a limit. The public MRI cohorts with repeated follow-up scans that I have not already used (for example Burdenko-GBM-Progression on The Cancer Imaging Archive, 180 patients) are under the NIH controlled-access policy. The request must be filed in dbGaP by a faculty-level researcher, and trainees may work on the data only under that person's oversight.

Would you consider sponsoring this request and overseeing my use of the data? In practice this means:

- you (with your institution's signing official) file the dbGaP request for accession phs004225;
- the data stay on storage you approve;
- I write a dated analysis protocol before any image is read, and you see it first;
- you review any result before it is shown to anyone.

I can send a two-page summary of the project, including the parts that did not work. So far one prediction test failed its pre-set rule and one imaging analysis is exploratory only. I would rather show you that openly.

I know this is a real request on your time. A short call or a "no" would both help me.

Thank you.

[name]
[school]
[contact]

---

## After sending

- Record the date sent and to whom in `START_HERE.md`.
- Do not download or open any shared data until a dated protocol for that dataset exists.
- Check with a parent or teacher before agreeing to any data-use terms.
