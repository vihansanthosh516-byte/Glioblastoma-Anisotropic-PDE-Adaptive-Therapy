# Track B ACRIN 6677: Stage 0 structure-only record (2026-11-12)

Protocol: `paper/trackB_acrin6677_protocol_2026-11-12.md` (frozen in commit 85378cb before download). This record reports structure, headers, and presence counts only. **No measurement value, progression date, survival time, or treatment arm was printed, summarized, or correlated.**

## File

- `data/acrin6677/ACRIN-DSC-MR-Brain-Full-Clinical.zip`, 2,303,087 bytes, SHA-256 `584a17f37017e5a00962a1677cc4640a6073f79ec1ae905f639e0867d9866157`, downloaded 2026-11-12 from the TCIA record (CC BY 4.0 as stated there). Not committed (`data/acrin6677/` is ignored).
- Contents: two folders set by the publisher, "Training-75" and "Test-25", each with form CSVs, `ACRIN_6677 Data Dictionary.xlsx`, and a form-description PDF. No licence file inside the zip.

## Findings

- The central-reader measurements are in form **S0 "CENTRAL READER QUANTITATIVE FORM"**: `s0e1` time point (Baseline, Week 2, Week 8, Week 16, ... Week 96, Off-study, Interval scan), `S0e2d` days from base date to MRI, `s0e4deid` reader ID, `s0e23` TOTAL TUMOR AREA (components labelled CM2), `s0e33` 3D TUMOR VOLUME, `s0e43` T2 MRI TUMOR VOLUME.
- S0 rows: 1,086 (856 + 230). Patients: 123. Distinct reader IDs: 46.
- Patients with at least three measured visits (non-missing value and date): 3D volume **70** (103 with any); 2D total area 83 (121 with any); T2 volume 90 (123 with any).
- Visits with two or more readers: 3D 394 of 415; 2D 430 of 490; T2 530 of 533.
- The form-description PDF does not contain the S0 form (text search for "tumor volume", "tumor area", "quantitative" returned nothing). The dictionary label for `s0e33` gives no unit.

## Gates

| Gate | Result | Evidence |
|---|---|---|
| S1 per-visit central-reader measurements with patient ID and time | **Pass** | Form S0 fields above |
| S2 at least 25 patients with three measured visits | **Pass** | 70 patients for 3D volume |
| S3 written definition of the measurement and its unit | **Not passed** | `s0e33` "3D TUMOR VOLUME" has no unit in the dictionary; the S0 form is absent from the PDF. Under the protocol, unknown is not a pass |
| S4 licence permits analysis and reporting | Pass | CC BY 4.0 on the TCIA record; no licence file inside the zip |

## Consequence

By the frozen rule, the analysis does not start. The gap is a documentation gap about a unit, not a data gap. Two honest paths:

1. A dated amendment, written before any value is read, that takes the unit from the trial's own publication (Boxerman et al. 2013, Neuro Oncol 15:945-954) and fixes the constant `c` from it. The analysis would then be labelled "pre-specified with one amendment made before any measurement value was read". The amendment must quote the paper's unit statement.
2. Stop and write a no-go.

Also noted for any amendment (structure facts, not results): the publisher's Training-75 / Test-25 split exists; the protocol uses all patients as one cohort and that stays. Time points include "Week 2" and "Interval scan" rows; the protocol's rule (first three measured visits by date) already covers them.
