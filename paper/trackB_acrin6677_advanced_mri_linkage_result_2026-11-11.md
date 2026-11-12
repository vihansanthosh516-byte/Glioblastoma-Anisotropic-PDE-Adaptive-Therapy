# ACRIN 6677 sRCBV linkage feasibility result

**Outcome: stop under the frozen rule.** The metadata audit found 15 eligible patients; the predeclared minimum was 25. No rCBV magnitude was evaluated, no model was fitted, and no patient score was produced.

## Frozen audit

- Protocol: [advanced-MRI increment protocol](trackB_acrin6677_advanced_mri_protocol_2026-11-11.md), frozen before table-row audit.
- Code: `src/track_b_acrin6677_advanced_mri_linkage.py`, committed as `deea04a`.
- Audit time: 2026-11-11 23:57 America/Chicago.
- Inputs: S0 and SC clinical tables from the Training-75 and Test-25 folders.
- Inspected: patient key, nominal visit, scan day, reader key, QC codes, and blank/nonblank status for the predeclared S0 area and SC sRCBV fields.
- Not evaluated: numeric area or rCBV magnitudes, treatment arm, survival, progression, or images.

The standard CSV parser decoded each row into temporary text fields before the script retained its allow-listed columns. It did not convert, display, summarize, or save non-allow-listed numeric fields. This is a semantic metadata-only audit, not a byte-level guarantee that other columns never entered a temporary parser buffer.

## Counts

| Check | Count |
|---|---:|
| S0 table rows | 1,086 |
| S0 rows excluded as uninterpretable | 29 |
| Patients with three valid S0 area visits | 85 |
| SC table rows | 192 |
| Patients with exact patient/visit/day SC linkage | 33 |
| Patients with an adequate SC row at that linkage | 15 |
| Patients with an adequate, nonempty sRCBV field | **15** |
| Frozen minimum | **25** |

All 15 eligible patients had nonempty presence flags in each of the three sRCBV lesion slots, and one qualifying reader row each. These are presence counts only; the values were not opened for analysis.

The audit stopped because 15 is below 25. The project did not relax the threshold, widen the date match, substitute nRCBV, or switch to DCE or spectroscopy. The advanced-MRI analysis is not run on this cohort under this protocol.

## Reproducibility record

- Aggregate output: [`output/track_b_acrin6677_advanced_mri_linkage_20261111/linkage.json`](../output/track_b_acrin6677_advanced_mri_linkage_20261111/linkage.json)
- Output SHA-256: `191A44CD1F9AB8A1FF566292600E2CC6AA45874640DBBC85F094B3CB437EA79A`
- Protocol SHA-256: `88AF51C555FAD32EEFFE6557ED17BECFF82E3E2AB91E7991A7ABAFD827C9DE64`
- Audit-code SHA-256: `B365498B4962998004188AF79B4980EFFA47FE1B99553C444FA24E62B3FE37A1`
- Synthetic tests: 7 passed in a scratch copy before the audit.

The JSON stores aggregate counts and a SHA-256 of the eligible patient-key set. It stores no patient IDs, reader IDs, measurement values, images, or patient-level rows.

## Interpretation

This is a feasibility no-go, not evidence that perfusion lacks predictive value. The exact-date, exact-visit, adequate-quality sRCBV subset is too small for the frozen model gate. Prior ACRIN 6677 studies already examined DSC rCBV changes against progression-free and overall survival; they report a small optional DSC cohort. The proposed next-area endpoint differs, but this result cannot support a new predictive claim. See [Schmainda et al. (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4490871/) and [Boxerman et al. (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9909012/).

**Next useful direction:** resume metadata-only screening for a different, independent longitudinal glioblastoma MRI cohort with a documented outcome measure and at least 25 eligible patients. Do not reopen ACRIN values or weaken this audit's rules to create a larger sample.
