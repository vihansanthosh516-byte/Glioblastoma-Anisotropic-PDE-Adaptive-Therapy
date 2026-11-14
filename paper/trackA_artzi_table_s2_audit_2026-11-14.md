# Track A Artzi supplementary metadata audit

Source: local publisher supplement cd-25-2302_table_s2_suppst2.xlsx, Sheet 1, header in row 1 and patient records in rows 2–41. The workbook has 40 patient metadata rows. It includes patient ID, age, sex, MGMT status, survival fields, primary tumor location, extent of primary surgery, and active TMZ at progression. Primary tumor location is a brain location, not a specimen region label.

The table has no specimen ID, surgery or time-point key per specimen, PIT or core category, pathology assignment, or malignant-cell annotation join. Therefore the number of same-patient, same-surgery PIT/core pairs is undetermined. Patient rows cannot be counted as eligible pairs. No expression data or patient outcomes were analyzed.

Harwood confirmed by email shared by the user that the older Zenodo version 2 annotations do not join to pathology and pointed to PIT in the new study. The next metadata source is Supplementary Table S1 or a specimen manifest. A useful manifest must link coded patient, resection/time point, specimen or field ID, independent pathology region and cell annotation key. Verify primary IDH-wildtype status and malignant-state coverage before test access. [Artzi et al. and supplementary listings](https://doi.org/10.1158/2159-8290.CD-25-2302).

## Supplementary Table S1 follow-up

Source: user-supplied cd-25-2302_table_s1_suppst1.pdf, page 1. Table S1 is an aggregate comparison of the Artzi study and two published cohorts. It gives a study total of 40 patients and contains no patient-by-specimen rows, PIT or core labels, pathology assignments, or same-surgery keys. The eligible pair count remains undetermined. A coded specimen or field-of-view manifest is needed. No test expression was opened.
