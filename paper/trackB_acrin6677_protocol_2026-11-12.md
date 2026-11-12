# Track B: ACRIN 6677 tabular forecast protocol (frozen 2026-11-12, before the data file is opened)

**Status: frozen before any ACRIN 6677 file was downloaded or opened.** Written from the public TCIA record only (https://www.cancerimagingarchive.net/collection/acrin-dsc-mr-brain/, read 2026-11-12). The commit that adds this file is the freeze record. Any later change is a dated amendment and makes the affected analysis exploratory.

> **Historical-status note (2026-11-12):** This initial protocol specified 3D tumor volume. Its volume-unit gate was not resolved; the analysis was amended to 2D area, then completed with a failed result. See `paper/trackB_acrin6677_amendment1_2026-11-11.md` and `paper/trackB_acrin6677_result_2026-11-11.md`. The later sRCBV feasibility audit also stopped at 15 of the required 25 patients; see `paper/trackB_acrin6677_advanced_mri_linkage_result_2026-11-11.md`.

This is a new Track B analysis on a cohort this project has never used. It does not change the Duke result (exploratory, spent), the LUMIERE no-go, Track A's failed A5 rule, or the Track C stop and its 14-day target.

## 1. Data

- Source: `ACRIN-DSC-MR-Brain-Full-Clinical.zip` (2.2 MB, CC BY 4.0) from the TCIA record above. Trial: RTOG 0625 / ACRIN 6677, randomized phase II, bevacizumab with irinotecan or temozolomide, recurrent glioblastoma. Record states 123 subjects and that two central neuroradiologists serially measured 2D and 3D enhancement on post-contrast T1 and the volume of FLAIR hyperintensity.
- Images are under NIH controlled access and are **not** used. No voxel is read.
- Stored outside Git under `data/acrin6677/`; the zip's SHA-256 is recorded at download.

## 2. Stage 0: structure-only inspection (no outcome values)

Allowed: list archive members; read data dictionaries, PDFs of blank forms, and column headers; count rows; count, per patient, how many visits have a non-missing central-reader measurement (missing or present only); record the visit-time variable's unit.

Not allowed in Stage 0: printing, plotting, summarizing, or correlating any measurement value, progression date, survival time, or treatment arm.

Stage 0 gates (all must pass, unknown is not a pass):

- S1. The open zip contains per-visit central-reader tumor measurements (3D enhancing volume, or 2D product if 3D is absent) linked to a patient ID and a visit time.
- S2. At least **25 patients** (the project's frozen minimum) have at least **three** measured visits (two history, one target) with the target visit after the second history visit.
- S3. A written data dictionary or form defines the measurement and its unit.
- S4. The licence text in or beside the file permits analysis and reporting (CC BY 4.0 as stated on the record).

If any gate fails: stop, write a no-go record, run nothing else on this cohort.

Stage 0 also fixes, in a dated addendum **before any value is read**: the exact variable names, which reader's value is used (rule below), and the visit labels. It may not change sections 3 to 6.

## 3. Question and estimand

In recurrent glioblastoma patients treated with bevacizumab-based therapy, does a two-parameter growth-with-treatment-effect curve, fitted to each patient's own first two measured visits plus the baseline rule below, forecast the tumor measurement at the patient's third measured visit with lower error than simple extrapolation?

- Unit: patient. One forecast per patient (first eligible window only: measured visits 1, 2 as history, measured visit 3 as target).
- Measurement: 3D enhancing volume if S1 finds it, otherwise 2D enhancing product. Fixed by availability in Stage 0, not by results.
- Reader rule: the mean of the two central readers when both exist; otherwise the single available reader.
- Scale: `y = log(measurement + c)`, with `c` = the smallest positive measurement unit stated in the dictionary (fixed in Stage 0 from the dictionary, not from the data distribution).

## 4. Methods compared (all fixed now)

All use only the patient's own two history visits `(t1, y1)`, `(t2, y2)` and the target time `t3`. No parameter is learned across patients, so there is no cross-patient leakage and no tuning.

- **B1 last value:** predict `y2`.
- **B2 log-linear:** extend the line through the two history points to `t3`.
- **B3 damped log-linear:** `y2 + 0.5 * slope * (t3 - t2)`.
- **Candidate M (growth with decaying treatment effect):** `dy/dt = g - k * exp(-t / tau)`, with `tau` fixed at 56 days (one 8-week assessment interval; fixed here, not fitted). `g` and `k` are solved exactly from the two history increments under the constraint `g >= 0`, `k >= 0`, using baseline `t = 0` at the first measured visit; when the exact solution violates a constraint, the violating parameter is set to 0 and the other is refit by least squares. This is the same structural form as Track C's tumor dynamics (net growth minus a treatment kill term) and is the only link to Track C: a check of that form on real trajectories, not a validation of the simulator.

Rationale stated honestly: bevacizumab often gives an early drop in enhancement followed by regrowth; M can represent that and the baselines cannot. M can also be worse, because two points barely determine two parameters.

## 5. Endpoint, test, and decision rule

- Per-patient error: `|predicted y3 - observed y3|`.
- Primary comparison: M versus the **best of B1, B2, B3**, where "best" is the baseline with the lowest median error on this cohort (this choice favours the baselines).
- Test: two-sided exact Wilcoxon signed-rank test on per-patient paired error differences, alpha 0.05.
- **Success** (both required): median paired difference favours M; p < 0.05.
- **Failure:** the best baseline is significantly better than M, or the median favours a baseline.
- **Inconclusive:** anything else.
- Reported regardless: n, medians with bootstrap 95% intervals (10,000 patient resamples, seed 20261112), the win count, and all three baseline results.

Secondary, descriptive only (no success claim): results by treatment arm; results using the other measurement type; FLAIR volume; later windows.

## 6. Synthetic tests before the patient run

A test file must pass in a scratch copy before any measurement value is read: (a) M recovers `g`, `k` on noise-free curves generated from M; (b) on pure log-linear data B2 beats M or ties; (c) constraint handling returns finite predictions for flat, rising, and falling histories; (d) the scoring code gives the known answer on a hand-built three-patient table.

## 7. Claim boundary

- Supports, if successful: on one fresh trial cohort, a simple treatment-effect growth form forecast the next central-reader measurement better than naive extrapolation.
- Does not support: drug resistance, biology-MRI integration (no expression data), clinical benefit, validation of the Track C simulator or its policies, or generalization beyond recurrent GBM on bevacizumab.
- Known weaknesses: enhancement under bevacizumab partly reflects vessel permeability, not tumor burden (pseudo-response); the trial authors already analysed early progression versus survival (Boxerman et al. 2013, Neuro Oncol 15:945-954), so the novelty of this forecast comparison is **not yet checked** and must be checked by a prior-art search before any novelty claim; a single window per patient; reader variability.

## 8. Not decided here

Nothing. Variable names and the constant `c` are fixed by Stage 0 from documentation only.
