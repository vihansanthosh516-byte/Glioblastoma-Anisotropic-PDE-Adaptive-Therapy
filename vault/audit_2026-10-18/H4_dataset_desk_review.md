# H-4 step 3: desk review of candidate locked-test datasets against Amendment A13.4

Date: 2026-10-18. Status: factual desk check from public description pages. No dataset was downloaded, no image or label was opened, and no dataset was chosen. Under A13.4 the choice is yours, on eligibility and logistics only, declared in an amendment before any label is read.

## Rules followed (plan v2.1, P3)
- Allowed: dataset description pages, licence terms, cohort and annotation documentation.
- Excluded: any page, table or section reporting model performance, benchmark results or leaderboards for any candidate.
- An 'unknown' cell is not inferred. An unknown on any A13.4 criterion blocks selection until you resolve it with the data provider.

## A13.4 criteria
C1. Pre-operative MRI with a usable core segmentation AND a recurrence MRI with enhancing / all-recurrence segmentation (labels 1-3, 3 = enhancing, 4 = cavity, or a documented mapping).
C2. At least 40 patients expected evaluable under A13.5 (time1 enhancing recurrence and standard-plan coverage < 0.999).
C3. Public or authorised access with a research-use licence.
C4. Not previously scored by this project in any form.

## Result
| Criterion | UCSD-PTGBM (TCIA) | BraTS 2025 Task 1 (adult glioma) | Duke longitudinal (Zenodo 20548528) |
|---|---|---|---|
| C1 pre-op + recurrence with enhancing label | UNKNOWN. Described as post-treatment (post-operative) imaging by design; pre-operative scans are not mentioned; annotation scheme not stated in the pages seen | UNKNOWN. A review describes Task 1 as combining pre- and post-treatment glioma data, with ET and a resection-cavity label for post-treatment scans; whether each patient has a pre-op scan AND a later recurrence scan is not stated | YES on the record: pre-operative and recurrence scans are identified by suffix per patient; labels 1 necrotic core, 2 non-enhancing T2/FLAIR, 3 enhancing, 4 resection cavity. The segmentations are automatic (nnU-Net), not expert. Whether that counts as 'usable' is your call |
| C2 at least 40 evaluable | UNKNOWN. 178 patients and 243 timepoints stated; 40 patients held out for the 2025 BraTS challenge; evaluable subset cannot be known without labels | UNKNOWN. A participant paper's data section gives 2,877 cases in total (unverified against the official page); number of paired pre-op + recurrence patients not stated | UNKNOWN. 61 patients and 503 sessions stated; evaluable subset cannot be known without labels |
| C3 access and licence | YES. CC BY 4.0 via TCIA | UNKNOWN. Licence for the glioma post-treatment data not stated in the pages seen; access route not confirmed | UNKNOWN. Files are restricted and, per a third-party index entry summarising the Zenodo record, released after a request with affiliation and purpose; the licence text was not seen |
| C4 never scored by this project | YES by the project's own record (H-2: TUM, MU, LUMIERE, RHUH; H-3: UCSF-PDGM / ALPTDG). Overlap with BraTS 2025 holdout not checked | UNKNOWN. The UCSD-PTGBM page says 40 of its patients are in BraTS 2025; whether BraTS 2025 also contains cohorts this project already used (for example RHUH or UCSF) is not stated and must be checked | YES by the project's own record (not among the sets used) |

## Reading
- No candidate satisfies all four criteria on public metadata. Every candidate has at least one unknown, which under the stated rule blocks selection.
- Duke is the only one with an explicit yes on C1 (pre-op and recurrence scans with the ET label scheme), but its labels are automatic, its access is by application, and its licence text was not seen.
- UCSD-PTGBM has a clear licence but its description points to post-treatment imaging only; C1 needs confirmation from the TCIA collection page or the authors.
- BraTS 2025 has the most unknowns, and an overlap risk with sets this project already used.
- The A13.4 choice rule ('earliest stable access wins' among qualifying candidates) cannot be applied until the unknowns are resolved.

## What would resolve the unknowns (all by you or the data providers, not by downloading)
1. Duke: submit the Zenodo access request and read the rights text; ask the authors whether recurrence labels are manual anywhere.
2. UCSD-PTGBM: read the TCIA collection page's data description for pre-operative scans and label scheme.
3. BraTS 2025: read the official challenge data page for licence, per-patient pairing and the list of contributing cohorts.
4. For any set: C2 can only be confirmed once labels are read, which A13.4 forbids before the amendment. The amendment should therefore state the expected evaluable count and how it was estimated from metadata.

## Exclusion log
Search results also returned pages that report model performance (for example an IDH foundation-model benchmark and several BraTS challenge participant papers on segmentation accuracy). None was read for results or used. The single number taken from a participant paper (2,877 cases) is a cohort size from its data section and is marked unverified.

Sources (description pages): UCSD-PTGBM Scientific Data / PMC12886825 and TCIA page; Zenodo record 20548528; BraTS dataset review (PMC11945730) and BraTS 2024/2025 data descriptions. Not independently verified beyond these pages.
