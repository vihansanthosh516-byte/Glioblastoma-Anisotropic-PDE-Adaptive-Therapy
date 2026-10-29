# Track B PDE external cohort screen (2026-10-28)

**Status:** metadata-only source review. No access request or contact was made. No restricted files, patient-level outcomes, images, masks, or scores were accessed.

## Strongest current lead: Duke longitudinal GBM MRI

The [Zenodo record](https://zenodo.org/records/20548528), *A Longitudinal Glioblastoma MRI Dataset with Anatomical Landmark Pairs for DIR Validation*, describes 61 patients with 3–18 sessions each (503 total), skull-stripped 1-mm isotropic atlas-aligned NIfTI, and T1, T1-c, T2, FLAIR, ADC, and DWI at each session. Relative scan and treatment dates are available. Every session has an nnU-Net mask with necrotic core, nonenhancing T2/FLAIR abnormality, enhancing tumor, and resection cavity labels. The record was published 2026-06-12 and modified 2026-09-23; its data files are restricted and total 119.3 GB.

The record offers CC BY-NC-ND 4.0 or a custom Duke license. A request asks for academic affiliation, if any, and a research-purpose description; access is not guaranteed. Before any request, confirm that the selected license and derivative-data rules fit the school, science fair, and intended outputs. No request was submitted.

## Fit and limitations

- Multi-visit masks make this a plausible lead for a **new, separately amended patient-level PDE growth-forecast evaluation**. A candidate design uses two historical masks to forecast the next eligible mask. The visit schedule is irregular, so freeze whether the estimand is next-observed-visit conditional on the realized gap or a prespecified fixed-day window supported by the available dates. Do not describe a realized-gap analysis as a scheduled-at-origin forecast.
- The published record says full cohort-selection details will appear in a future paper. Its recurrence category can be assigned by tumor-volume change or clinical progression, so a recurrence-conditioned analysis can select on tumor change. The masks are automated nnU-Net outputs, not independent expert ground truth. Treat a mask-based endpoint as released-label change unless validated against an independent reference.
- The 61 patients are not a power guarantee. Before outcome access, audit patient overlap and prior exposure against MU, UCSD-PTGBM, CFB-GBM, and other project cohorts; estimate precision under the final patient-level estimand; and freeze patient-grouped splits, baseline comparators, parameter-identifiability and calibration checks, and technical-failure rules.
- UCSD-PTGBM is a separate, partially exposed annotation-change side study with a negative exploratory result; it is not the PDE core. Duke would require a new protocol/amendment and cannot retroactively change that result.
- CFB-GBM is a separate annotation-transition candidate with unresolved review lineage and registration QC. A CFB result cannot complete the PDE parameter-estimation and forecasting program.

## Decision

Duke is the best metadata-level lead found for a future multi-visit Track B PDE study, but it is **not ready for scoring**. First resolve permitted use and access, then audit overlap/exposure, registration, and target-label limitations. If these gates fail, retain the current exploratory evidence; do not substitute another cohort after inspecting outcomes.

## Sources

- Primary dataset record and license: [Zenodo 10.5281/zenodo.20548528](https://zenodo.org/records/20548528).
- Current Track B boundaries and CFB preflight: `START_HERE.md` and `paper/track_b_cfb_scoring_protocol_draft_2026-10-27.md`.
- Historical cohort and access audit: `README_STATUS.md` section 5a.1 and ledger entry 2bm.
