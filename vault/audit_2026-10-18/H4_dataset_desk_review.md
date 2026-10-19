> **HISTORICAL BACKGROUND (corrected 2026-10-19):** Dated old-GBM audit record, not current instructions or authorization to run analyses. Current goal/status: START_HERE.md. Simulation settings are assumptions, not measured treatment effects. Original approval, chronology and process claims are author-reported unless supported by retained logs.

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
| C2 at least 40 evaluable | UNKNOWN. 178 patients and 243 timepoints stated; Version 3 partitions are 136 subjects/184 studies and 42 subjects/59 studies in BraTS-GLI 2024 test data; earlier 40/BraTS 2025 description was wrong; evaluable subset cannot be known without labels | UNKNOWN. The former 2,877 participant-paper count is omitted without an official scoped cohort count; number of paired pre-op + recurrence patients not stated | UNKNOWN. 61 patients and 503 sessions stated; evaluable subset cannot be known without labels |
| C3 access and licence | YES. CC BY 4.0 via TCIA | UNKNOWN. Licence for the glioma post-treatment data not stated in the pages seen; access route not confirmed | CONDITIONAL. Official Zenodo Rights lists CC BY-NC-ND 4.0 or a custom Duke agreement; access requires a request. Approval and permitted use remain unresolved |
| C4 never scored by this project | YES by the project's own record (H-2: TUM, MU, LUMIERE, RHUH; H-3: UCSF-PDGM / ALPTDG). Overlap with BraTS 2025 holdout not checked | UNKNOWN. The current TCIA record describes 42 subjects in BraTS-GLI 2024 test data; whether BraTS 2025 also contains cohorts this project already used (for example RHUH or UCSF) is not stated and must be checked | YES by the project's own record (not among the sets used) |

## Reading
- No candidate satisfies all four criteria on public metadata. Every candidate has at least one unknown, which under the stated rule blocks selection.
- Duke is the only one with an explicit yes on C1 (pre-op and recurrence scans with the ET label scheme), but its labels are automatic, its access is by application, and its stated rights/access conditions still require approval and permitted-use review.
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

Source correction checked 2026-10-19 (metadata only, ledger 2bm): official Duke record https://zenodo.org/records/20548528 documents pre-op/recurrence suffixes, automatic masks and labels 1-4; official TCIA Version 3 https://www.cancerimagingarchive.net/collection/ucsd-ptgbm/ documents 178 subjects/243 timepoints and CC BY 4.0. These totals do not establish 40 evaluable cases. BraTS source lead https://zenodo.org/records/15094823 does not establish pairing/license eligibility. Before any future selection, an amendment must define a metadata-only estimate procedure for evaluability and its uncertainty; no label inspection is authorized here.
