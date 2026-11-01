# Track B fixed-growth baseline: draft protocol wording (2026-10-31)

Status: **candidate wording, resolved choices (revised 2026-10-31); not frozen until the support gate passes**. No patient data was scored, no restricted file was accessed, no model was run. This file adds a baseline definition and a search update. It does not authorize scoring. Sources checked: `paper/trackABC_execution_plan_2026-10-31.md` (B2.4–B2.7, commit `5ab4e86`) and `paper/research_plan_complete_project_appendix.md` §3.2–3.5.

## 1. Candidate protocol wording (supersedes the first draft)

Code: `src/fk_fixed_growth.py`, wired into `src/fk_mu_patient_runner.py`. Synthetic tests: `tests/test_fk_fixed_growth.py` (17 passed); with the MU runner, patient-CV and precision-gate suites, 65 passed (scratch copy, Python 3.14, 2026-10-31). No patient mask was read.

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

## 3a. Runner integration (2026-10-31)

- `src/fk_fixed_growth.py:add_fixed_growth_arm` wraps the patient-CV callbacks. `src/fk_mu_patient_runner.py` now calls it, so every MU development forecast also carries the fixed-growth mask, `dice_fixed_growth`, and `paired_delta_dice_hybrid_minus_fixed_growth`. The result adds `fixed_growth_paired_difference_bootstrap` (seed 20261101).
- The radius is chosen inside the outer-fold nuisance fit, which `run_outer_patient_cv` calls with outer-training patients only. Training patients with an empty visit-2 core are left out of selection and counted.
- The forecast is clipped to the loader's origin-only domain (`origin_domain_by_patient`), not to target support. Scoring then applies the same rule as the other arms: a forecast outside validated target support is an error.
- The runner still refuses to run unless the loader records both source-backed gates (registration and target support) as passed. This change does not touch those gates. No gate is passed.

## 5. Frozen fold-training recipes (candidate; freeze date 2026-10-31)

Both recipes see outer-training patients only. Neither sees future tumor size. Neither is tuned on outer-test patients.

**U-Net (required mask baseline).**
- Inputs: visit-1 core, visit-2 core and origin-only domain on the 2 mm grid (three channels). These are the same inputs the PDE arm uses. An MRI-intensity variant is not part of the required comparison.
- Model: 3D U-Net, 4 resolution levels, 16 base channels doubled per level, instance norm, LeakyReLU, one output logit.
- Training: one visit-2 → visit-3 pair per outer-training patient. Augmentation: random flips on each axis and random 90° rotations in the axial plane. Loss: soft Dice + binary cross-entropy, equal weight. Adam, learning rate 1e-3, batch 1, 300 epochs, last epoch kept. No early stopping and no inner tuning: there are too few patients to tune without overfitting the selection.
- Seeds: 20261102 + outer fold index, three seeds per fold. Average the three probability maps, threshold at 0.5, clip to origin-only domain.
- Empty visit-2 core: copy-last (empty), as for the other arms.
- Failure: a non-finite loss or an all-empty training output means "U-Net unavailable"; the novelty test is then **incomplete**, not passed.
- Known limit: about 14 training patients per fold is small for a U-Net. Report this next to the result.

**Training-only count model (count metrics only).**
- Count: number of 26-connected components of visit-3 core minus the visit-2 core grown by 2 mm, keeping components of at least 2 voxels (16 mm³). This is a label-change count, not a biological lesion count.
- Model: Poisson GLM, log link, predictor log(visit-2 core volume in mm³), offset log(actual interval in days). Fit on outer-training patients only. If the fit fails or does not converge, use the training mean rate per day times the interval.
- Scores: Poisson deviance, count MAE, and Brier score for P(count ≥ 1). Compare against the training mean-rate model. **Never scored with Dice and never turned into a mask.**
- Status: coded in `src/fk_learned_baselines.py` and wired into the MU runner (2026-10-31).

**U-Net fixed settings (added 2026-10-31, before coding; `src/fk_learned_baselines.py:UNET_RECIPE`).**
- Preprocessing: three binary float channels (visit-1 core, visit-2 core, origin-only domain) on the 2 mm grid. No intensity data, no normalization. Zero-pad each axis to a multiple of 8, and pad the two in-plane axes (0, 1) to the same size so axial 90° rotations are defined. Crop the output back to the original shape.
- Architecture: encoder widths 16, 32, 64, 128 (four levels). Each level has two 3×3×3 convolutions with padding 1, each followed by affine instance norm and LeakyReLU (slope 0.01). Down-sampling: 2×2×2 max pool. Up-sampling: 2×2×2 transposed convolution, then concatenation with the skip and two convolutions as above. Output: 1×1×1 convolution to one logit.
- Loss: soft Dice (smoothing 1.0) + binary cross-entropy with logits, equal weight, over the whole padded volume. The target is the visit-3 core clipped to the origin-only domain.
- Optimizer: Adam, learning rate 1e-3, betas (0.9, 0.999), weight decay 0. Batch 1. 300 epochs. Each epoch visits every outer-training pair once in an order shuffled by the seeded generator. Augmentation per step: flip each axis with p = 0.5, then rotate the axial plane by k × 90°, k uniform in {0, 1, 2, 3}. Same transform on inputs and target. Keep the last epoch.
- Seeds: 20261102 + outer fold, + 1000, + 2000 (three models). Deterministic torch algorithms on; CPU float32.
- Prediction: sigmoid of each model, mean of the three probability maps, threshold ≥ 0.5, clip to origin-only domain. Empty visit-2 core: copy-last (empty).
- Unavailable: non-finite loss, or every training-pair prediction empty after training. Then the U-Net arm is "unavailable" for that fold, the rows record it, and the novelty test is incomplete.
- Training patients with an empty visit-2 core are left out of training and counted.

**Proposed, not adopted:** a history-centred 64³ U-Net crop and U-Net reuse across prior-support runs (draft amendment A-2026-10-31-U in `paper/track_b_runtime_pilot_2026-10-31.md`). The code supports the crop only through `UNetRecipe.crop_voxels`, which is `None` in the frozen recipe.

**Count model fixed settings (`src/fk_learned_baselines.py`).**
- New label piece: a 26-connected component of (visit-3 core) \ (visit-2 core dilated by 2 mm), with at least 2 voxels. Both masks clipped to the origin-only domain first.
- Model: Poisson GLM log μ = b0 + b1 log(visit-2 core volume, mm³) + log(interval days). Bug fix 2026-10-31: if the training total count is 0, the mean rate uses 0.5 in the numerator. This keeps a predicted count above 0, because a predicted count of 0 gives infinite deviance for any new piece. Exact rule: mean_rate = Σy/Σt if Σy > 0, else 0.5/Σt, over outer-training patients only. It is training-only smoothing, frozen before any patient outcome; do not tune it against test counts. Report how many folds used it. Fit by Newton–IRLS, at most 50 iterations, tolerance 1e-8 on the coefficient change. Fallback to the training mean rate per day if fewer than 3 usable training patients, a singular step, non-convergence, or non-finite coefficients.
- Blind forecast: the forecast step outputs a rate per day from history only. The scorer multiplies it by the actual interval. The forecast never sees the target day.
- Empty visit-2 core: use the training mean rate.
- Scores per patient: Poisson deviance, absolute count error, Brier score of P(count ≥ 1) = 1 − exp(−μ). The training mean-rate model is scored the same way as the reference. No mask and no Dice.

**Distance-shell spatial rule.** Dropped on 2026-10-31. It was never a required baseline: the roadmap listed it as optional "only if frozen before scoring", and it was not frozen. The roadmap now says dropped. No required baseline was removed.

## Remaining blockers

- Source-backed valid support S and MRI-only longitudinal registration QC (MU and Duke). Scoring is blocked until both exist.
- Cavity meaning from source label documents.
- Duke access and school-use approval.
- U-Net and count model: coded and wired into `src/fk_mu_patient_runner.py` (`src/fk_learned_baselines.py`); 13 synthetic tests passed, and 78 with the fixed-growth, MU runner, patient-CV and precision-gate suites (scratch copy, Python 3.14). Runner tests use a tiny U-Net recipe for speed; the runner itself always passes the frozen recipe. The frozen U-Net (300 epochs × 3 seeds × 5 folds, CPU) has **no timing pilot yet**; runtime is unknown. Until the U-Net runs on an eligible held-out cohort, the novelty test is incomplete.
- Fixed growth is wired into the MU development runner only. The Duke runner does not exist.

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
