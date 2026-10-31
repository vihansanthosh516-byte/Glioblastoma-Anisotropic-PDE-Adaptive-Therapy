# Track B fixed-growth baseline: draft protocol wording (2026-10-31)

Status: **candidate wording, resolved choices (revised 2026-10-31); not frozen until the support gate passes**. No patient data was scored, no restricted file was accessed, no model was run. This file adds a baseline definition and a search update. It does not authorize scoring. Sources checked: `paper/trackABC_execution_plan_2026-10-31.md` (B2.4–B2.7, commit `5ab4e86`) and `paper/research_plan_complete_project_appendix.md` §3.2–3.5.

## 1. Candidate protocol wording (supersedes the first draft)

Code: `src/fk_fixed_growth.py`. Synthetic tests: `tests/test_fk_fixed_growth.py` (12 passed, scratch copy, Python 3.14, 2026-10-31). No patient mask was read.

**B-FG. Fixed-growth comparator.** For a held-out patient with last history tumor-core mask C on the forecast grid, approved valid-support mask S, and voxel spacing h (mm), the forecast is

F_r = { v ∈ S : d_h(v, C) ≤ r }

d_h is the exact Euclidean distance (mm) from voxel v to the nearest voxel of C, computed with spacing h (`scipy.ndimage.distance_transform_edt(~C, sampling=h)`). C has d_h = 0. r = 0 gives C ∩ S, which is copy-last clipped to S.

1. **Horizon.** One nominal day-90 window. Eligible pairs have an actual history-to-target interval of 76–104 days. Report each actual interval. The forecast is conditional on the pair falling in that window. It is not exact-day-90 validation. One window means one radius per outer fold.
2. **Grid and spacing.** Use the existing protocol grid: 2 mm isotropic, nearest-neighbour resampling after the registration gate passes (`paper/track_b_pde_external_protocol_draft_2026-10-28.md` §3). All mask comparators use this grid.
3. **Radius grid.** R = {0, 2, 4, …, 20} mm: steps equal one 2 mm voxel; finer steps would not change the mask. The cap is 20 mm, ten voxels. It is a search bound, not a biological claim. If the chosen radius equals 20 mm, record "boundary hit". Do not widen the grid after scoring starts.
4. **Radius selection.** For each outer fold, pass only that fold's outer-training patients to `select_radius`. Score each r by plain patient-macro Dice of F_r against the training patients' full target cores. Pick the best mean. Ties within 1e-9 go to the smaller r. **No inner folds.** The comparator has one free choice and no fitted parameter, so inner folds would choose the same r and add nothing. Outer-test targets never enter selection; the forecast for a test patient is built before its target is read (`run_outer_fixed_growth`).
5. **Nonnegative radius.** The grid has no negative radius. Fixed growth never erodes. If r = 0 is chosen, report "fixed growth = copy-last" for that fold. It is not a separate win or loss.
6. **Empty history.** An empty last core is ineligible: the code raises an error. This matches the existing runner, which sends empty history to the persistence fallback, and appendix §4.3, which requires a nonempty input core. Count such patients in cohort flow.
7. **Empty target and Dice.** Use `dice_with_empty_rule` from `src/fk_mu_scoring.py`: both empty = 1, one empty = 0. Score against the full target core. Target core outside S stays in the target and counts as a miss for every method; report its voxel count per patient.
8. **Valid support S.** S is an input. This module does not define it. S must be a source-backed valid-acquisition-support mask mapped to the forecast-origin grid, with MRI-only registration QC that does not use GTVs. The earlier draft's rule ("union of nonzero MRI voxels and labels") is withdrawn: no source establishes it. **Until a source-backed S exists, patient scoring stays blocked.** The same S clips every mask comparator, including copy-last.
9. **Cavity.** The cavity is not removed from S by default. Remove it only if the source label documentation says target tumor-core labels cannot occur in a prior cavity. Decide this from source documents only. Do not inspect target-mask voxels to decide it. If the meaning stays unclear, scoring stays blocked.

**B-CMP. Comparators (same held-out patients, same origins, same targets, same window, same S, same grid):**

| Comparator | Role | Score |
|---|---|---|
| Copy-last (persistence), clipped to S | required mask baseline | plain Dice |
| Fixed growth (B-FG) | required mask baseline | plain Dice |
| PDE-only | required mask baseline | plain Dice |
| U-Net trained on outer-training patients only | required mask baseline | plain Dice |
| Training-only count model | count model, not a mask | count metrics only |
| Volume-trend (matched-volume signed distance), geometric train-rate | existing comparators (§3.3) | plain Dice, as now |
| Component ablations | existing comparators (§3.3) | plain Dice, as now |

The distance-shell mask rule is **dropped** for now (see §3).

**Pass rule.** The candidate passes only if its patient-macro plain Dice beats **every required mask baseline** under the frozen margin and interval rule. No model and no scorer receives future tumor size. Volume-matched comparators that draw sizes from training-only models stay allowed; oracle volume matching does not.

## 2. Appendix reconciliation (done in this commit)

`paper/research_plan_complete_project_appendix.md` §3.3 and §3.5 now list the four required mask baselines, the count model on count metrics, the dropped distance-shell rule, and the "beat every required mask baseline" pass rule. The volume-trend, geometric-rate and ablation comparisons are kept. This amendment is dated 2026-10-31, before any Track B target mask is opened.

## 3. Choices resolved

| # | Choice | Decision | Reason |
|---|---|---|---|
| 1 | Horizon | One 76–104-day window, nominal day 90; actual interval reported; conditional estimand | Matches the existing Duke and MU design; avoids tiny bins |
| 2 | Cavity in S | Not removed by default; remove only if source label documents say core cannot occur in a prior cavity; never decided from target voxels | Avoids a rule that the data have not established |
| 3 | Radius grid | 0–20 mm in 2 mm steps; boundary hit reported | One-voxel steps on the 2 mm grid; cap is a search bound, not a radiation margin |
| 4 | Erosion | Not allowed | Removes a free choice; shrinkage is covered by the volume-trend comparators |
| 5 | Distance-shell rule | Dropped | Its fitting and failure rules are not frozen; keeping it would add an untested required baseline |
| 6 | Grid | 2 mm isotropic, nearest-neighbour | Same grid as the existing protocol |

## Remaining blockers

- Source-backed valid support S and MRI-only longitudinal registration QC (MU and Duke). Scoring is blocked until both exist.
- Cavity meaning from source label documents.
- Duke access and school-use approval.
- `src/fk_fixed_growth.py` is not yet wired into `src/fk_mu_patient_runner.py`. Wiring it in is a code task. It does not authorize scoring.
- U-Net and count model: no frozen training recipe yet.

## 4. Search update (continues `paper/track_b_prior_art_followup_2026-10-31.md`)

Evidence levels: FT full text, AB abstract, SN snippet, API metadata.

**PubMed (E-utilities, 2026-10-31; title/abstract fields):**

| Query | Hits | Relevant finding |
|---|---|---|
| glioma/GBM AND "new lesion(s)" AND count/Poisson/negative binomial/number AND predict/forecast | 2 | Neither is a count forecaster (PMID 35419519 volumetric thresholds; 33952802 circulating tumor cells). |
| glioma/GBM AND reaction-diffusion/Fisher-Kolmogorov/proliferation-invasion AND deep learning/neural network/physics-informed/hybrid AND MRI | 8 | PMID 40389199 (Zhong et al., J Theor Biol 2025, HybrSyn): growth model + deformable registration to synthesize MRI; scored by SSIM 80.93% on BraTS19 + one hospital (AB). Not a mask forecast against later scans. Others are inverse/infiltration work already known (LMI, PINN infiltration). |
| glioma/GBM AND growth AND Dice AND persistence/"no change"/"last observation"/"previous scan"/dilation | 0 | No PubMed paper with a persistence-Dice comparison. |
| glioma/GBM AND growth model/reaction-diffusion AND segmentation/registration AND uncertainty/perturbation/robustness/sensitivity | 3 | PMID 37058375 (Liang et al., IEEE TMI 2023): Bayesian tissue heterogeneity, 4 rats, Dice 0.89, accuracy depends on number of earlier scans (AB). No patient perturbation study. |

**Google Scholar: not done.** Scholar has no public API and blocks automated queries. This search gap stays open.

**Citation sweep (Semantic Scholar API, 2026-10-31):**
- PREDICT-GBM (2509.13360): 2 citing records, OASA (2609.12435, known) and an MRI-compression paper. Nothing new.
- GliODIL (2312.05063): 33 citing records. New names to read: "A Differentiable Brain Tumor Mass Effect Solver" (2610.06097), "Physics-Grounded Disentangled Flow Modeling for Brain Disease Progression Trajectory" (2606.28630), "An MRI-informed poromechanical model for organ-scale prediction of glioma growth" (2605.05013), "Brain Tumor Growth Inversion via Differentiable Neural Surrogates" (no ID).
- TaDiff (2309.05406): 37 citing records. New names: "Uncertainty Gated Residual Learning with Patient Level Bagging for Longitudinal Glioma Progression Prediction" (2026), "3D-GlioPREDICT" (2606.05113), "Brain-WM: Brain Glioblastoma World Model" (2603.07562), "MUMINS" (2609.17169).
- SoC-DT: 0 citing records. LUMIERE, MU-Glioma-Post and UCSF-ALPTDG dataset searches hit API rate limits (HTTP 429). **Not done.**

**Full-text and abstract checks:**
- **arXiv 2607.13877 (FT, first 100k characters):** the residual has a learnable scale initialized to 0.1 and a hard cutoff |residual| < 0.1 → 0 during rollout. Weight decay 1e-5. No explicit penalty toward the PDE. Anatomy from UPENN-GBM, but all longitudinal tumor evolution is simulated. Dice one-step: RD 0.641 ± 0.189, dense hybrid 0.920 ± 0.057. No code URL found.
- **arXiv 2606.28630 (AB):** Wang and Liu. PDE-regularized flow matching (diffusion-reaction-advection) on three public longitudinal brain-disease datasets, not named on the abstract page. No metrics on the abstract page. Code github.com/jhuldr/PDF (not checked).
- **"Uncertainty Gated Residual Learning…" (SN only, UNVERIFIED):** the only text found was on a ResearchGate page with wrong metadata. It describes a radiomics progression classifier with a residual correction applied to low-confidence predictions, and patient-grouped CV. This is a classification task, not a spatial forecast. But "gating a residual by uncertainty" is close in idea to the stability-shrunk residual. Authors, venue and data are unconfirmed. **Must be read before any claim about the residual.**

**What this changes.** Nothing found here reports plain Dice against copy-last by horizon. Nothing found is a glioma count forecaster. The residual-gating idea now has at least one possible neighbour (unverified). The gap is still preliminary; no "first" claim.

**Still open:** Google Scholar; dataset citation sweeps (rate-limited); full texts of 2606.28630, 2610.06097, 2605.05013, 3D-GlioPREDICT, Brain-WM and the uncertainty-gated residual paper; SoC-DT code once released.
