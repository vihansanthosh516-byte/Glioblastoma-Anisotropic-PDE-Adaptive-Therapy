# Track B code and data audit (2026-11-14) — post hoc, exploratory

Scope: the next-scan core-volume forecast. Question, target and score unchanged.
Scripts and outputs: `output/track_b_code_audit_20261114/` (`audit1`, `audit2`, `fix_and_loco4`, `loco4_all_baselines`, each with `_stdout.txt`).
All four cohorts (MU, UCSD, UCSF, LUMIERE) are spent for this question. Nothing in this file is a test result.

## 1. What was checked and found correct

| Check | Result | File |
|---|---|---|
| MU label meaning | 1 = non-enhancing tumor core, 2 = FLAIR hyperintensity, 3 = enhancing, 4 = resection cavity. Core = 1 + 3 is correct. Source: https://www.nature.com/articles/s41597-025-06011-7 | — |
| MU voxel counts in the pair table vs raw masks | 0 mismatches in 387 pairs | `audit1_stdout.txt` |
| MU scan days and intervals vs the clinical table | 0 mismatches | `audit1_stdout.txt` |
| MU pair = next scan in time | 0 wrong | `audit1_stdout.txt` |
| MU voxel size and grid | 1 mm, 240x240x155, all 594 masks; no unexpected label | `audit1_stdout.txt` |
| UCSF mask voxel count vs the dataset's own volume table | ratio 1.0 (5% to 95%: 1.0 to 1.0) | `audit2_stdout.txt` |
| Bootstrap interval code | matches normal-theory interval on toy data | `audit2_stdout.txt` |
| Per-patient averaging | correct on a toy case | `audit2_stdout.txt` |

No bug was found that changes a reported number by a large amount.

## 2. Problems found

| # | Problem | Size | Fix | Status |
|---|---|---|---|---|
| 1 | UCSD: a session is dropped when "days to initial surgery" is "Unknown". 12 of 49 multi-scan patients lose their pairs. | 46 pairs / 33 patients become 55 / 42 | Take the interval from any event fixed in time (initial surgery, death, last follow-up, first radiation, first chemo). Anchors agree: median difference 0 days; 4 intervals differ by more than 7 days, the median is used. | Fixed in `fix_and_loco4.py` (new code; the locked script is unchanged) |
| 2 | LUMIERE: 39 scans have a mask and no rating row. One first scan labelled `week-000-1` has no rating, so the locked rule did not remove it as pre-operative. | 1 scan | Treat an unrated `week-000-1` scan as pre-operative. | Fixed in `fix_and_loco4.py`; pair count stays 400 |
| 3 | MU: one patient (PatientID_0081) has a second scan day of -965. | 1 pair, already excluded by the interval rule | None needed | Recorded |
| 4 | UCSF: 14 eligible rows have impossible day values; treatment windows were assumed. | Affects the treatment-phase input only | Range check | Recorded in `paper/trackB_problem_review_2026-11-14.md` |
| 5 | MU: the 32 worst pairs hold 36.5% of the error. | — | See section 3 | Not a code problem |

## 3. What the worst MU pairs are (not a code problem)

32 pairs with |log change| > 3 (`audit1_stdout.txt`, section D):
- 20 falls. In 9 the resection cavity grew by more than 1,000 voxels (a later surgery is likely). In 8 the FLAIR label grew by more than half of the lost core (enhancement went away while FLAIR increased; this is the pattern of bevacizumab response or of relabelling). In 5 the whole tumor fell by more than 1 log.
- 12 rises. The input core was tiny (1 to 219 voxels in 9 of 12) with a large FLAIR region. Enhancement then appeared.
- These are real events in the masks. They are not visible in the input size alone.

Other cohorts (`audit2_stdout.txt`):
- UCSD: 16 pairs with radiotherapy between the scans hold 46% of the error.
- LUMIERE: 99 pairs that start at a post-operative scan hold 33% of the error. 22 pairs that end in a post-operative scan hold 5%.

## 4. After the fixes: leave-one-cohort-out over four cohorts (exploratory)

Each cohort is held out in turn; every method is fitted on the other three. 442 patients. `loco4_all_baselines_stdout.txt`.

Pooled patient MAE: VIG 1.121, persistence 1.148, population rate 1.154, Gompertz (absolute loss) 1.189, Gompertz 1.217, EB slope 1.252, gradient boosting (absolute loss) 1.307, gradient boosting 1.435, two-point 1.680.

VIG minus each comparator, pooled, 95% CI: persistence -0.027 (-0.039, -0.014); population rate -0.032 (-0.046, -0.019); Gompertz -0.095 (-0.139, -0.052); Gompertz absolute loss -0.067 (-0.107, -0.029); EB slope -0.131 (-0.165, -0.098); gradient boosting -0.313 (-0.396, -0.232); gradient boosting absolute loss -0.186 (-0.254, -0.118); two-point -0.558 (-0.882, -0.342). Every upper limit is below 0.

VIG minus persistence by held-out cohort: MU -0.036 (-0.061, -0.012); UCSD -0.019 (-0.063, +0.023); UCSF -0.024 (-0.039, -0.008); LUMIERE -0.023 (-0.058, +0.011). VIG is the only method that is not worse than persistence in any cohort (worst case -0.019; Gompertz worst case +0.152, population rate +0.039).

### Weakest points (read before quoting)

- **This is not a confirmed win.** The pooled rule was chosen after the locked LUMIERE run failed its own rule (`output/track_b_vig_lumiere_20261114/result.json`). All four cohorts are spent.
- The form of VIG (interval curve, rise/fall history terms) was chosen with the outcomes of MU, UCSD and UCSF in view. Leave-one-cohort-out refits the numbers but not the design. Only the LUMIERE fold is free of design selection, and alone it does not pass (CI crosses 0).
- The gain is small: 0.027 of 1.148, about 2.4%.
- In single cohorts VIG does not beat the best baseline of that cohort: LUMIERE Gompertz is +0.020 better; MU Gompertz (absolute loss) is level (+0.003).
- The interval signal is clinical behaviour, not biology (see the protocol).
- Gradient boosting here uses only the inputs present in all four cohorts.
- The pooled interval treats patients from four cohorts as one sample.

## 5. Prior-art check (partial)

Two web searches, result snippets only. "Informative observation times" is a known idea in biostatistics (for example PSA testing frequency in prostate cancer follow-up, arXiv 1508.07511). No paper was found that uses the scan interval as a predictor of glioma volume at the next MRI. This is NOT a full prior-art search. Novelty is unverified.

## 6. What a confirmed win needs

One untouched cohort, the same frozen VIG code, and the same "beat all comparators" rule fixed before the run. Using the pooled effect (-0.027) and the per-patient SD seen on LUMIERE (0.154, `output/track_b_vig_lumiere_20261114/result.json`), about 7.85 x 0.154^2 / 0.027^2 = 255 patients give 80% power against persistence (my calculation; it assumes the same effect and spread).
