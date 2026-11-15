# Track B: screen for an untouched test cohort (2026-11-15)

Purpose: find data that can give the VIG hybrid one fair, locked test. Question, target and score are unchanged.
Counts: `output/track_b_cohort_screen_20261115/screen_stdout.txt` (script `screen.py` in the same folder).

**What was read:** metadata only (scan days, file lists, dataset pages, clinical columns known before the outcome scan).
**What was NOT read:** no mask, no image, no tumor volume at an outcome scan, no radiomics table. No candidate cohort is spent by this screen.

## 1. What a test cohort needs

1. Two or more post-treatment scans per patient, with the days between them.
2. Core volume (necrosis label + enhancing label) at both scans, or images from which one frozen tool can make those labels.
3. Access without a senior sponsor.
4. Outcome volumes not yet seen in this project.
5. About 255 patients for 80% power against persistence (my calculation; effect 0.027 from `loco4_all_baselines_stdout.txt`, SD 0.154 from `output/track_b_vig_lumiere_20261114/result.json`; it assumes both carry over).

## 2. Result of the screen

| Cohort | Access | Patients with an eligible pair (upper limit) | Labels | Outcome seen before? | Verdict |
|---|---|---|---|---|---|
| **WAW** (Zenodo 10.5281/zenodo.22205631) | Open, CC BY 4.0 | 28 (48 pairs; 25 are grade 4) | Same four labels as MU: 1 necrotic core, 2 FLAIR, 3 enhancing, 4 resection bed. FL-PoST/FeTS then radiologist-refined. | No. Only four metadata tables are on disk. | **Fits the target exactly. Too small alone.** |
| **CFB-GBM v4** (TCIA) | Open, CC BY 4.0. Already on disk. | 158 (155 t0-to-t1 pairs, 107 t1-to-t2 pairs) | One binary GTV label only (radiotherapy target). No core labels. | No. Zero GTV voxels were read (`START_HERE.md`, CFB section). | **Largest. Needs core labels from one frozen tool.** |
| RHUH-GBM (TCIA) | On disk | About 11 (estimate from the clinical table) | Expert sub-region labels | **Yes.** Recurrence masks were read in the H-2/H-3 test (`src/123_predictgbm_test_score.py`). | Too small and exposed. Not used. |
| Burdenko-GBM-Progression (TCIA) | **NIH Controlled Data Access** (TCIA page, checked 2026-11-15) | 180 patients in total | GTV | No | Blocked. Needs a senior sponsor. |
| Duke (Zenodo 20548528) | Custom licence; school approval pending (`START_HERE.md`) | 61 patients in total | Automated labels | **Yes.** Target masks were read on 2026-11-04. | Used. Not a fresh test. |
| BraTS 2024 post-treatment (Synapse) | Registration | Unknown | Same four labels | Four of seven sites are our spent cohorts (MU, UCSF, UCSD, Duke) | No scan dates were found in the challenge paper (arXiv 2405.18368, passages read). Not usable without intervals. Unverified for the full data release. |
| ACRIN 6677 | — | — | — | **Yes** (area forecast failed, 2026-11-11) | Spent. |
| PREDICT-GBM, UPENN-GBM | On disk | — | — | — | Input scan is pre-operative. Does not fit the question. |
| MU, UCSD, UCSF, LUMIERE | — | — | — | **Yes** | Spent. |

Other sources checked and rejected from their abstracts: Yale-Brain-Mets-Longitudinal (metastases), IMAGO, UTSW glioma, OpenNeuro ds007045 (no sign of serial post-treatment scans; unverified), private cohorts in method papers (Oslo, Heidelberg).

## 3. Plain result

- **No single open cohort is large enough.** The largest clean pool is CFB + WAW: 186 patients at most.
- Power against persistence with 186 patients is about 0.67 (`screen_stdout.txt`). With WAW alone it is about 0.15. The "beat all eight comparators" rule is harder than this single comparison, so real power is lower.
- These counts are upper limits. The "input core > 0" rule can only be applied after the masks exist.

## 4. Problems to solve before any test

| # | Problem | Fix | Weakest point of the fix |
|---|---|---|---|
| 1 | CFB has no core labels. Using its GTV would change the target. | Train one segmentation model on spent cohorts (MU, UCSF, UCSD masks). Inputs: T1Gd + FLAIR (527 CFB scans have both; 430 also have pre-contrast T1; only 144 have T2). Freeze it. Then label CFB. | The target becomes our own tool's output. Tool error adds noise. Standard four-sequence tools cannot run because T2 is missing. |
| 2 | Mask source differs between CFB (our tool) and WAW (released masks). | Primary analysis on released WAW masks and tool-made CFB masks, cohort as a stratum. Sensitivity: run the same tool on WAW images. | Two mask sources in one pooled test. |
| 3 | CFB t0 is before radiotherapy. Radiotherapy lies between t0 and t1. | Keep these pairs (UCSD has the same kind). Report t0-to-t1 and t1-to-t2 separately as secondary. | Large treatment effect in the first pair; persistence may be poor there for every method. |
| 4 | CFB scan times are in whole weeks. | Use weeks x 7. | Coarser interval than in other cohorts. |
| 5 | CFB intervals follow a plan (4 and 6 months). | None needed. The spread is real: 35 to 245 days and 14 to 210 days (`screen_stdout.txt`). | The interval may carry less of the "doctor is worried" signal than in routine care. This is the part of VIG that did not show on LUMIERE. |
| 6 | Power is below 80%. | State it before the run. A failed test is then "not confirmed", not "disproved". | A pass at 67% power is still a pass; a fail is weak evidence. |
| 7 | One try only. | Lock code, parameters, rule and patient list before any CFB label or WAW mask is opened. | — |

## 5. Options

**A. Pooled locked test on CFB + WAW (recommended).**
Steps: (1) write and commit the protocol; (2) train and freeze the segmentation model on spent cohorts (GPU, Kaggle); (3) download `3_2_postop_masks.zip` (2.6 MB) from Zenodo; (4) label CFB; (5) run frozen VIG once.
Cost: several days, mostly the segmentation model. Chance of a pass: below 67% by the power estimate, and VIG failed its one locked test before.

**B. Small locked test on WAW only.**
Steps: protocol, download 2.6 MB, run once.
Cost: hours. Power about 0.15. It would spend the only cohort that matches the target exactly, for almost no evidence. Not recommended as a confirmation test.

## 6. Weakest points of this screen

- One web search round and one paper-index round. A cohort can exist that was not found.
- BraTS 2024 date availability was judged from the challenge paper, not from the data files.
- RHUH count uses progression-free time as a stand-in for the scan interval.
- The power numbers assume the pooled exploratory effect holds. That effect was chosen post hoc (`paper/trackB_code_audit_2026-11-14.md`).
