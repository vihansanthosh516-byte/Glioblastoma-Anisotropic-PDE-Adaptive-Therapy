# Track B PDE external cohort screen (2026-10-28)

**Status:** metadata-only source review. Duke declined the first access request because planned-use and academic-affiliation details were missing. The user reports resubmitting on 2026-10-29 with the requested research-purpose and license-compliance information; access approval and school/fair permitted-use approval are pending. No restricted files, patient-level outcomes, images, masks, or scores were accessed.

## Strongest current lead: Duke longitudinal GBM MRI

The [Zenodo record](https://zenodo.org/records/20548528), *A Longitudinal Glioblastoma MRI Dataset with Anatomical Landmark Pairs for DIR Validation*, describes 61 patients with 3–18 sessions each (503 total), skull-stripped 1-mm isotropic atlas-aligned NIfTI, and T1, T1-c, T2, FLAIR, ADC, and DWI at each session. Relative scan and treatment dates are available. Every session has an nnU-Net mask with necrotic core, nonenhancing T2/FLAIR abnormality, enhancing tumor, and resection cavity labels. The record was published 2026-06-12 and modified 2026-09-23; its data files are restricted and total 119.3 GB.

The record offers CC BY-NC-ND 4.0 or a custom Duke license. A request asks for academic affiliation, if any, and a research-purpose description; access is not guaranteed. Duke declined the initial request for missing details, and the user reports resubmitting on 2026-10-29 with the requested information and license-compliance statement. School/fair acceptance of the selected license and derivative-data rules, and access approval, remain unconfirmed.

## Fit and limitations

- Multi-visit masks make this a plausible lead for a **new, separately amended patient-level PDE growth-forecast evaluation**. A candidate design uses two historical masks to forecast the next eligible mask. The visit schedule is irregular, so freeze whether the estimand is next-observed-visit conditional on the realized gap or a prespecified fixed-day window supported by the available dates. Do not describe a realized-gap analysis as a scheduled-at-origin forecast.
- The published record says full cohort-selection details will appear in a future paper. Its recurrence category can be assigned by tumor-volume change or clinical progression, so a recurrence-conditioned analysis can select on tumor change. The authors' [dataset guidance](https://github.com/deshanyang/GBM-DIR-QA) adds that separate sessions have not been registered to one another, and recommends rigid or affine registration before downstream use. It says scan dates were extracted automatically from exam/report dates and may be imperfect or off by a few days; its recurrence workflow can label a scan after >40% enhancing-volume growth from the postoperative nadir. These create registration, timing, and outcome-selection risks for a longitudinal forecast. The masks are automated nnU-Net outputs, not independent expert ground truth. Treat a mask-based endpoint as released-label change unless validated against an independent reference.
- The 61 patients are not a power guarantee. Before outcome access, audit patient overlap and prior exposure against MU, UCSD-PTGBM, CFB-GBM, and other project cohorts; estimate precision under the final patient-level estimand; and freeze patient-grouped splits, baseline comparators, parameter-identifiability and calibration checks, and technical-failure rules.
- UCSD-PTGBM is a separate, partially exposed annotation-change side study with a negative exploratory result; it is not the PDE core. Duke would require a new protocol/amendment and cannot retroactively change that result.
- CFB-GBM is a separate annotation-transition candidate with unresolved review lineage and registration QC. A CFB result cannot complete the PDE parameter-estimation and forecasting program.

## PREDICT-GBM is not a fresh cohort

The new [PREDICT-GBM platform](https://doi.org/10.1038/s41746-026-03194-0) reports 243 patients across TUM-GBM (142), LUMIERE (61), and RHUH-GBM (40), and releases processed data and an evaluation pipeline. Those cohorts already appear in this project’s development or prior analyses: the local PREDICT-GBM work used the 142 TUM patients, and LUMIERE/RHUH were also examined. Treat the platform as useful for reproducibility and baseline review, not as untouched external validation. Its primary endpoint is geometric coverage of later recurrence for radiotherapy planning, which is related to but distinct from a PDE forecast of future tumor-core geometry. See the local exposure record in `README_STATUS.md` and the [PREDICT-GBM data release](https://huggingface.co/datasets/LZimmer/PREDICT-GBM).

## Decision

Duke is the best metadata-level lead found for a future multi-visit Track B PDE study, but it is **not ready for scoring or download**. The user reports that an access request is submitted; confirm school/fair acceptance of the chosen license and derivative rules while access is pending. If access is granted, audit overlap/exposure, session registration and landmark QC, date reliability, and target-label limitations before scoring. If these gates fail, retain the current exploratory evidence; do not substitute another cohort after inspecting outcomes.

## Sources

- Primary dataset record and license: [Zenodo 10.5281/zenodo.20548528](https://zenodo.org/records/20548528).
- Dataset-specific registration and date guidance: [GBM-DIR-QA README](https://github.com/deshanyang/GBM-DIR-QA).
- Current Track B boundaries and CFB preflight: `START_HERE.md` and `paper/track_b_cfb_scoring_protocol_draft_2026-10-27.md`.
- Historical cohort and access audit: `README_STATUS.md` section 5a.1 and ledger entry 2bm.
