# Data-request plan (2026-11-12)

**Decision (user, 2026-11-12): data-access requests are now allowed.** This replaces the project's earlier no-request rule (gate G5 in `paper/trackB_fresh_cohort_screen_2026-11-11.md`). No other gate, frozen rule, endpoint, or saved result changes. Track C confirmatory evaluation stays stopped and its 14-day target is unchanged (user decision 2026-11-10).

Nothing was downloaded, no request was sent, and no one was contacted in this pass. All facts below were read from the cited pages on 2026-11-11 or 2026-11-12 (metadata only).

## Why requests are needed

Open, no-request sources were screened and none gives a patient-level link among tumor biology, repeated MRI, and measured drug response:

| Source | Finding | Evidence |
|---|---|---|
| HGCC (Johansson 2020, Cell Reports) | Expression is public (GSE152160, 145 samples, HTA 2.0 arrays). No measured temozolomide (TMZ) response found: TMZ appears only in the introduction and reference list of the main text; MGMT is absent. Table S2 (drug list) not opened, so absence is likely, not proven. | https://www.cell.com/cell-reports/fulltext/S2211-1247(20)30878-0 ; GEO GSE152160 |
| Ntafoulis 2023 (Br J Cancer) | Measured TMZ response on 66 cultures (agent-reported, paper not opened by the main session); RNA-seq for 24 cultures (GEO GSE232173, GPL20301, verified). Paired count is at most 24; too small for a defensible held-out model. | GEO GSE232173 |
| Ivy GAP (TCIA) | 39 patients, 390 MRI studies (pre-surgery, post-surgery, follow-up). Region-labelled RNA-seq: 122 samples from 10 tumors only. Images under NIH Controlled Data Access Policy. No follow-up tumor labels listed. BraTS 2021 is built partly on this collection. | https://www.cancerimagingarchive.net/collection/ivygap/ |

## Eligibility barrier (most important fact)

TCIA controlled-access collections are requested through dbGaP (phs004225 for most TCIA face-containing datasets). The TCIA policy page states that requestors "must be at a level equivalent to a tenure-track professor, or senior scientist", and that trainees "are not permitted to submit data access requests in dbGaP" but "may be part of projects using such data when the projects are overseen by an eligible user". Source: https://www.cancerimagingarchive.net/nih-controlled-data-access-policy/ (read 2026-11-12).

**Consequence: the author cannot file the dbGaP request alone. A sponsoring senior researcher (professor or senior scientist with an institutional signing official) must file it and oversee the project.** Whether the author has such a sponsor is not recorded anywhere in the project.

## Ranked request targets

| Rank | Target | What it would give | Who must request | Limits |
|---|---|---|---|---|
| 1 | Burdenko-GBM-Progression (TCIA) | 180 patients, 1 to 8 follow-ups each; MGMT and IDH status; treatment-response status (response, progression, pseudoprogression); RT target volumes (GTV/CTV/PTV). Fresh Track B cohort. | Sponsor via dbGaP | No FLAIR or enhancing-tumor labels, only RT volumes, so the frozen Track B endpoint cannot be reused unchanged; a new dated protocol is required. No expression data. MGMT is a proxy, not measured drug response. 131 GB. |
| 2 | Ivy GAP images (TCIA) plus Ivy GAP clinical database (registration) | Patient-level link between region-labelled RNA-seq and repeated MRI. | Sponsor via dbGaP; separate registration at ivygap.org | At most 10 tumors with region RNA: pilot only, cannot meet the 25-patient minimum. RNA already used in Track A (exposed). No tumor labels. Per-patient visit counts unverified. 140 GB. |
| 3 | Moniz Garcia et al. 2026 (Neuro-Oncology, Mayo Clinic Florida) | 6 patients: MRI-guided multi-sampling, bulk RNA-seq and exome, and measured drug tests incl. TMZ (CellTiter-Glo) in patient-derived 3D models. The only source found with biology, MRI, and measured drug response in the same patients. | Email to corresponding authors; statement: "Sequencing data are available upon request ... to comply with institutional policy" | n = 6: case series only. Longitudinal MRI not stated. Institutional agreement likely needed. https://pmc.ncbi.nlm.nih.gov/articles/PMC13550674/ |
| 4 | Ntafoulis 2023 (Erasmus MC) | Per-culture TMZ values and the tissue RNA-seq (n = 56 "on reasonable request", agent-reported) to enlarge the paired set beyond 24. | Email to corresponding author | Still small; in-vitro response only. |

## What each could and could not support

- Burdenko: a fresh, pre-registered Track B progression analysis. Not drug resistance; not biology-MRI integration.
- Ivy GAP: an exploratory, pre-registered pilot of expression-versus-MRI-change association in the same patients (n <= 10). Not confirmatory.
- Moniz Garcia: descriptive patient-level illustration of all three data types (n = 6). Not a statistical test.
- None of these changes Track A's failed A5 rule, Track B's Duke status (exploratory, spent), or Track C's stop.
- No positive result is guaranteed by any request.

## Before any data arrives

1. Identify a sponsor (blocking for ranks 1 and 2).
2. Draft request texts; the user and sponsor send them. Claude does not send requests.
3. Write a dated protocol per dataset before any image or value is read, with synthetic tests first.
4. Save every accepted data-use agreement in `paper/` (the Duke terms gap, U1, shows why).

## Not verified

Per-patient study counts for Ivy GAP; contents of the Ivy GAP clinical database; whether Burdenko follow-up scans include FLAIR for every visit; Ntafoulis data-availability sentence in full; dbGaP processing time.
