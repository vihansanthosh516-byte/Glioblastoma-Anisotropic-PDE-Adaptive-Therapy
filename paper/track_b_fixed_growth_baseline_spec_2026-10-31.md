# Track B fixed-growth baseline: draft protocol wording (2026-10-31)

Status: **draft, not frozen**. No patient data was scored, no restricted file was accessed, no model was run. This file adds a baseline definition and a search update. It does not authorize scoring. Sources checked: `paper/trackABC_execution_plan_2026-10-31.md` (B2.4–B2.7, commit `5ab4e86`) and `paper/research_plan_complete_project_appendix.md` §3.2–3.5.

## 1. Proposed protocol wording

**B-FG. Fixed-growth comparator.** For a held-out patient with last history tumor-core mask C (origin space, physical voxel spacing h), valid support S, and declared horizon bin b, the forecast is

F_r(b) = { v ∈ S : d_h(v, C) ≤ r(b) }

where d_h is the Euclidean distance in millimetres from voxel centre v to the nearest voxel of C, computed with the voxel spacing (C itself has d_h = 0). r(b) is one radius per horizon bin, chosen from the fixed grid R = {0, 1, 2, …, 15} mm.

1. **Inputs allowed.** The last history mask, the origin brain support, the origin cavity mask, voxel spacing, and the declared horizon bin. Nothing from the target scan enters the forecast: no target mask, no target volume, no target-derived support.
2. **Valid support S.** S = B_0 \ Q_0, where B_0 is the origin operational brain support (union of nonzero voxels in the origin MRI modalities and any nonzero origin label) and Q_0 is the origin resection cavity. The same S clips every mask comparator. Report target core voxels outside S as a count per patient. They count as misses for every method.
3. **Radius selection.** Inside each outer fold, use outer-training patients only. Split them into K = 5 patient-grouped inner folds with a seed fixed before scoring (fewer patients than 10 → leave-one-patient-out). For bin b and each r ∈ R, compute plain patient-macro Dice on the inner-validation patients whose eligible pair falls in bin b. Pool across inner folds. Pick r(b) with the highest mean inner Dice. Apply r(b) once to the outer-test patients in bin b.
4. **Tie rule.** If two radii are within 1e-6 of the best mean inner Dice, pick the **smaller** radius. So a tie with copy-last resolves to r = 0.
5. **Shrinking and zero growth.** The grid has no negative radius. Fixed growth never erodes. r = 0 is allowed and makes F equal copy-last clipped to S. If r(b) = 0, report "fixed growth = copy-last in bin b"; do not count it as a separate win or loss.
6. **Small bins.** If bin b has fewer than N_min = 5 outer-training patients, use the radius chosen on all outer-training patients pooled across bins, and flag bin b as "pooled radius".
7. **Empty cases.** Patients with an empty last core are not eligible (appendix §4.3 requires a nonempty input core; adopt the same rule for §3). If the target core is empty and a forecast is empty, Dice = 1; if exactly one is empty, Dice = 0. Use the same convention for every method.
8. **Reporting.** Report r(b) for each outer fold and bin, plus the spread of radii chosen within inner folds. Report plain Dice per patient, patient-macro Dice per bin, and paired patient-bootstrap differences against the candidate.

**Note on inner folds.** Fixed growth has no fitted parameter besides r. So pooled inner-fold Dice equals Dice over all outer-training patients in bin b, and the inner split does not change r(b). Keep it anyway: it makes every comparator use the same split, and the per-inner-fold radii show how stable the choice is.

**B-CMP. Comparator roles (same held-out patients, same origins, same targets, same bins, same S):**

| Comparator | Role | Score |
|---|---|---|
| Copy-last (persistence) | required mask baseline | plain Dice by bin |
| Fixed growth (B-FG) | required mask baseline | plain Dice by bin |
| PDE-only | required mask baseline | plain Dice by bin |
| U-Net trained within outer-training folds | required mask baseline | plain Dice by bin |
| Training-only distance-shell rule | separate mask baseline, only if frozen before scoring; then required | plain Dice by bin |
| Training-only count model | count model only | count metrics only (deviance, count MAE, calibration of P(count ≥ 1)) |
| Volume-trend, geometric-rate | existing appendix §3.3 comparators | as in appendix |
| Component ablations | existing appendix §3.3 comparators | as in appendix |

The candidate passes only if it beats every required mask baseline under the frozen margin rule (execution plan B2.6–B2.7).

## 2. Check against the current appendix

| Item | Appendix now | Gap |
|---|---|---|
| Mask comparators | §3.3: persistence, volume-trend, geometric-rate, ablations | No U-Net, no fixed growth, no distance-shell rule. **Amendment needed.** |
| Count model | not in §3 | Add as count-only secondary. **Amendment needed** (the execution plan already says this; the appendix does not). |
| Horizon bins | §3.2: realized interval, conditional estimand; Duke candidate is one 76–104-day window; §4.1 uses ≤90, 91–180, >180 days for the side study | No bins are declared for §3. **Must be frozen.** Bins use scan dates only, so the estimand stays conditional on the realized interval. |
| Valid support | §4.2 defines B_0 and Q_0 for the side study only | Adopt for §3 or define a §3 support. **Amendment needed.** |
| Tuning on training patients | §3.2: tuned comparator settings use outer-training patients only | Consistent. Add the inner-fold rule and tie rule. |
| Empty-core and Dice convention | §4.3 requires nonempty input core; §3 silent | **Amendment needed.** |
| Primary endpoint | §3.3: patient-macro Dice of the tumor core | Consistent. Keep. |
| Pass rule | §3.5: beat the strongest prespecified simple baseline | Execution plan B2.7 is stricter: beat every required mask baseline. **Reconcile in the amendment.** |

**Amendment needed before scoring:** one dated Track B amendment that (a) adds U-Net, fixed growth and the optional distance-shell rule to §3.3, (b) adds the count model as count-only secondary, (c) declares §3 horizon bins and N_min, (d) declares S and the empty-case Dice rule, (e) sets the pass rule to "beat every required mask baseline", and (f) fixes seeds and inner K. Freeze it before any target mask is opened.

## 3. Unresolved choices (decide before freezing)

1. **Horizon bins.** Duke-style single window (76–104 days) or three bins (≤90, 91–180, >180 days). A single window avoids tiny bins but tests one horizon only.
2. **Cavity in S.** Excluding Q_0 assumes the tumor core does not grow into the cavity. If the source labels show core inside the later cavity, keep Q_0 in S and report it.
3. **Radius grid.** 1 mm steps to 15 mm match the standard 15 mm margin used by PREDICT-GBM. A finer grid or a larger maximum would need a reason before scoring.
4. **Erosion.** Allowing negative radii would let fixed growth model shrinkage. It also adds a free choice. The draft does not allow it.
5. **Distance-shell rule.** Keep it or drop it. If kept, freeze its shell width, threshold grid and minimum addition count, and it becomes required.
6. **Resampling.** All comparators must run at one frozen grid (native or 1 mm isotropic). Distances in mm follow the grid spacing.

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
