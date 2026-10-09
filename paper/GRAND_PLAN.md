# Grand plan: answer every blocker (2026-10-07)

Written 2026-10-07. Source list: the 27 blockers + "not done" list given by the author today.
Every number cites a file in `output/` or `vault/PAPER_FINDINGS_LEDGER.md` (ledger IDs like "2w").
Literature tags: **[S]** = found by web search today (summary read, full paper not read; check before citing). **[VERIFY]** = detail not confirmed.
Compute: Docker Desktop on this PC (CPU, WSL2) + Kaggle GPU. This PC has an Intel Arc GPU, not NVIDIA, so CUDA work goes to Kaggle.

---

## 0. Read this first (the honest rule)

There are two kinds of blockers.

1. **Bugs and gaps.** We can fix these. Examples: 2 mm grid, clamp, 2,000 bootstraps, grid edge, W27 count.
2. **Findings.** These are true answers about the data. We cannot "fix" them. We can only ask a better question.
   Examples: the PDE does not beat persistence (ledger 2w), it does not beat a growth rule (2w), anisotropy adds nothing (2w, 2t), shape is close to noise (2s).

Tuning the PDE until it beats persistence on the same MU data is **not allowed**. A reviewer will find it. It is p-hacking.

So the plan has three parts:
- **Part A. Fix every bug.** Then re-run once, with the final settings, and report what comes out.
- **Part B. Turn the findings into a result.** "How much of next-scan GBM change can anyone forecast?" with a *measured* noise floor. This is a real, defensible result.
- **Part C. Test the PDE where theory says it should help.** That is not next-scan shape. It is **where the tumour comes back after surgery** (pre-op scan → recurrence). A public 3-centre benchmark exists for this (PREDICT-GBM, 253 patients). Its authors say anisotropic diffusion is **not modelled** yet. That is an open gap we can test.

"Revolutionary" here means: a pre-registered, public-benchmark answer to an open question. It may be positive or negative. Both are publishable if done right. Nothing in this plan promises a positive result.

---

## 1. The three headlines we aim for

| ID | Question | Data | Status |
|---|---|---|---|
| **H-1** (keep) | How forecastable is next-scan GBM change? Does any model beat persistence beyond the measured noise? | MU-Glioma-Post (105 GBM as analysed, ledger 2x), LUMIERE, RHUH-GBM | Results exist (2r, 2s, 2w). Needs measured noise floor + final run. |
| **H-2** (new primary, pre-register first) | Does a tract-anisotropic growth model cover more of the later recurrence than (a) the standard 15 mm margin and (b) an isotropic model, at the same target volume? | PREDICT-GBM: 253 patients, 3 centres, MIT licence, recurrence masks already in pre-op space | Not started. Must write Amendment 6 before any run. |
| **H-3** (new secondary) | Can we forecast **volume** change and next-scan progression with honest (calibrated) uncertainty? | MU (develop) → LUMIERE (frozen external test) | Not started. The oracle with true volume gives +0.128 Dice (2r), so volume is where the signal is. |

Why H-2 is the right place for the PDE:
- Reaction-diffusion models exist to estimate tumour cells **beyond the visible edge**. Recurrence after surgery comes from those cells. Next-scan shape after treatment is driven by surgery, radiotherapy and cavity change, which our PDE does not model.
- PREDICT-GBM (arXiv 2509.13360) [S]: mean recurrence coverage, aggregated cohort: standard 15 mm margin 77.34%, GliODIL 78.91%, their U-Net 79.37%, PINN-GBM 77.99%, LOTI 78.54%. All models are at equal volume to the standard plan. The paper states anisotropic diffusivity is not modelled and lists DTI as future work.
- A small lower-grade glioma study (14 patients) reported that a DTI-informed anisotropic model predicted recurrent shape better than isotropic (from a 2024 review) [S] [VERIFY original paper]. No large GBM test was found.

Risk, stated now: on MU, anisotropy changed Dice by +0.0002 (2w) and DTI alignment was null, p = 0.32 (2t). The PREDICT-GBM models beat the margin by only about 1.5 to 2 points. So H-2 may also be null. The pre-registration makes a null result count.

---

## 2. Compute setup (Docker + Kaggle GPU)

### 2.1 Docker Desktop (this PC)
**Installed and tested 2026-10-08** (Docker Desktop 4.94.0, WSL2 backend; `docker run hello-world` printed "Hello from Docker!"). Fix needed on the way: a leftover non-admin `C:\ProgramData\DockerDesktop` folder blocked the installer; deleted from an admin shell, then reinstalled.
Remaining step: test the project image with `make build && make test` (uses the existing `Dockerfile`) in P0.

What Docker is for here:
- Reproducible runs of the whole pipeline (`Dockerfile`, `docker-compose.yml`, `Makefile` already exist).
- Packaging **our** growth model in the PREDICT-GBM container format (input `/mlcube_io0`, output `/mlcube_io1`, see github.com/BrainLesion/PredictGBM). A CPU image is fine for our PDE.
- Limit: the PREDICT-GBM model images are built on NVIDIA CUDA. They will not use the Intel Arc GPU. **We do not need to run them**: the Hugging Face dataset already holds their prediction maps (`<model_id>_pred.nii.gz`).

### 2.2 Kaggle GPU
- Kaggle notebooks **cannot run Docker**. GPU code must be plain Python (PyTorch). Folder `kaggle_run/` already exists.
- Port the 3D solver (`run_improved_aniso.TensorFK`) to PyTorch so it runs on the GPU. Check it against the CPU solver on 5 cases before use (max voxel difference < 1e-4).
- Upload inputs as a **private** Kaggle dataset. Never upload patient data to a public dataset. Check each dataset licence allows this (PREDICT-GBM is MIT; TCIA collections have their own terms) [VERIFY MU-Glioma-Post terms].
- Weekly GPU quota and session length are limited [VERIFY current numbers on Kaggle]. Save checkpoints every patient. Write one output JSON per patient so a stopped session can resume.

### 2.3 Data to download (needs disk; 1.1 TB free on C: today)
| Data | Size | Source | Use |
|---|---|---|---|
| PREDICT-GBM processed | 25.5 GB | huggingface.co/datasets/LZimmer/PREDICT-GBM (MIT) | H-2 |
| SRI24 or other DTI tensor atlas in the same space as PREDICT-GBM | small | [VERIFY which atlas space PREDICT-GBM uses] | H-2 anisotropic tensor |
| LUMIERE imaging | 31 GB | already on disk, `data/external/lumiere/Imaging-v202211.zip` (unzip needed) | blockers 12, 18, H-3 |
| RHUH-GBM NIfTI | 2.9 GB | already on disk, `data/external/rhuh/` | noise floor (expert-corrected masks), H-2 overlap check |
| GBmap (single-cell, 240 patients) | large | Ruiz-Moreno et al., Neuro-Oncology 2025 [S] [VERIFY download portal] | Track A replication |
| Ivy GAP region RNA-seq | small | Allen Institute [VERIFY] | Track A region labels |

Command for PREDICT-GBM (run only after Amendment 6 is committed):
```
huggingface-cli download LZimmer/PREDICT-GBM --repo-type dataset --local-dir data/external/predict_gbm
```

---

## 3. Every blocker: type, fix, done-when

Type key: **Fix** = bug or gap, can be repaired. **Finding** = true result, report it; the action is a better test, not a tweak. **Accept** = cannot be fixed in this project; state it.

### A. Blocks a positive result

| # | Blocker | Type | What we do | Literature / data | Output | Done when |
|---|---|---|---|---|---|---|
| 1 | PDE does not beat persistence (GBM-only +0.0130, CI −0.0043 to +0.0302; 2w) | Finding | (a) Final run after all Part-A fixes (#14–19), 10,000 draws. (b) Report it as H-1. (c) Do **not** re-tune on MU. | — | `output/pde_manifest/results_final.json` | Final run done once; number in ledger, whatever it is |
| 2 | PDE does not beat growth rule (−0.0019, CI −0.0176 to +0.0141; 2w) | Finding | Add the **volume-matched shape test**: run PDE and geometric rule both scaled to the *true* next volume. Any gap left is pure shape skill. This tells us if the PDE knows anything about *where* the tumour grows, separate from *how much*. | Oracle idea from `output/baseline_ladder.json` | `output/shape_skill_matched.json` (new script 104) | Paired CI for PDE-shape minus geometric-shape at matched volume |
| 3 | Anisotropy adds nothing (+0.0002; 2w), DTI alignment p = 0.32 (2t) | Finding (on MU) | Keep as result. Re-test in the setting theory predicts: H-2 (pre-op → recurrence, longer horizon, infiltration beyond the edge). Use a tensor atlas; use patient DTI where it exists. Also try Painter–Hillen tensor build and free-water-corrected FA as pre-declared variants (max 2 variants, Holm). | PREDICT-GBM [S]; Painter & Hillen 2013 J Theor Biol [VERIFY]; free-water FA improved recurrence separation, AUC 0.9 vs 0.77 (PMC12696142) [S]; DTI review, 16 studies, 394 patients (PMC11626419) [S] | `output/predict_gbm/results.json` | H-2 run done under Amendment 6 |
| 4 | Shape nearly unforecastable (persistence Dice 0.277; 57% of pairs SNR < 1; 2s) | Finding | Make it H-1. Replace the assumed noise floor with a **measured** one (see #18). Report the forecastability ceiling per SNR bin. | BraTumIA test-retest: repeatability coefficient 31–46% for enhancing tumour, 20 patients scanned 2 days apart (AJNR 42(6):1080) [S]; inter-rater overlap poor after surgery (GCI 0.32) [S] | `output/forecastability_measured.json` | Ceiling recomputed with measured floor |
| 5 | Gain only on first pairs (+0.028), not later pairs (2w) | Finding + test | Pre-declared test of the "fitting effect" idea: (a) PDE vs geometric on first pairs is +0.012, CI −0.005 to +0.029 (2w), so even there the PDE adds no clear gain over growth. (b) Re-select cells on later pairs only, test on first pairs (swap). (c) Stratify by days since surgery and RT. | `oof_pairs.csv`, script 102 strata | `output/pde_manifest/first_vs_later.json` | Swap test reported |
| 6 | Signal is in volume (oracle +0.128; 2r); nothing predicts volume | Gap | H-3: volume forecast models. Ladder: persistence → population rate (geometric, 2r) → mixed-effects Gompertz → hierarchical Bayesian growth with treatment-phase covariates (days since RT, TMZ) → gradient boosting on pre-forecast features. Score: mean \|ln ratio\| error, CRPS, interval coverage. Patient-level conformal intervals. Secondary: RANO progression at next scan (LUMIERE AUC). | Conformal under shift needs weighting: coverage fell to about 80% vs 90% target without it (Lambert et al. MICCAI 2024, arXiv 2407.19938) [S]; patient-level conformal on sparse lesion series (arXiv 2609.21197, preprint) [S]; BraTS 2025 LUMIERE response task, one team AUC 0.81 (arXiv 2509.06511, preprint) [S]; earlier LUMIERE growth-vs-no-change result (2l) | `output/volume_forecast/results.json` | Dev on MU, frozen, then one LUMIERE run |
| 7 | Median gain tiny (+0.0095; W4) | Finding | Report full per-patient distribution, catastrophic-error rate (2v), and an **abstain rule**: forecast only when SNR > pre-set cut; otherwise say "change below detectable level". Score coverage vs accuracy. | Uncertainty-aware decision idea from conformal-band work (NeurIPS 2025, arXiv 2511.13911) [S] | `output/selective_forecast.json` | Accuracy-vs-coverage curve |

### B. Data problems

| # | Blocker | Type | What we do | Literature / data | Output | Done when |
|---|---|---|---|---|---|---|
| 8 | Small cohort (105 GBM) | Fix (partly) | Add cohorts, never double-count: MU 105 + RHUH 40 + LUMIERE 91 for H-1/H-3; PREDICT-GBM 253 for H-2. PREDICT-GBM includes 61 LUMIERE and 40 RHUH patients, so list overlaps by ID. Add a power simulation: minimum detectable Dice gain at our n. | PREDICT-GBM patient lists on GitHub [S] | `data/manifests/cohort_overlap.csv`, `output/power_sim.json` | Overlap table + MDE stated in paper |
| 9 | Informative follow-up: progressed 88.8% vs 47.8% (`output/cohort_tables.json`) | Fix (sensitivity) + Accept | (a) Inverse-probability-of-eligibility weights from baseline features (age, grade, MGMT, extent of resection). (b) Re-run primary with weights. (c) Wording: "among patients with usable scan pairs". IPW only corrects for *observed* causes; say so. | IPW + intensity weighting for informative dropout (arXiv 2510.19154) [S]; joint models when dropout depends on current tumour size [S] | `output/ipw_sensitivity.json` | Weighted and unweighted results side by side |
| 10 | 28 of 133 not GBM | Fixed by rule | Amendment 5 made GBM-only primary. Keep all-eligible as sensitivity. Nothing new. | — | — | Already done |
| 11 | 133 analysed vs 134 in manifest (W27, 2x) | Fix (document) | Amendment 6 states the pair-table rule. Regenerate `consort_mu.json` from code with both counts and the 3 patients named (0007, 0242, 0249). No hand edits. | — | `data/manifests/consort_mu.json` (regenerated) | One CONSORT figure, both numbers explained |
| 12 | LUMIERE masks machine-made; volumes noisy (2l) | Fix (measure) | Use the two automatic segmenters (HD-GLIO-AUTO, DeepBraTumIA) as two "raters". Their disagreement = measured noise. Check against RHUH **expert-corrected** masks. Label all LUMIERE results "automated masks". | LUMIERE paper: automatic tools need expert correction for response assessment [S]; RHUH-GBM expert-corrected masks (PMC10551826) [S] | `output/seg_noise_floor.json` | Floor in mm³ and ln units per volume bin |
| 13 | No prospective data | Accept | Best substitutes: (a) frozen config + git tag before any external run; (b) temporal split on MU (train on earlier-dated patients, test on later); (c) three external cohorts. Wording: "retrospective, in silico". | — | `configs/development.yaml`, tag `v1.0-development-frozen` | Tag exists before external run |

### C. Method and numerical problems

| # | Blocker | Type | What we do | Literature / data | Output | Done when |
|---|---|---|---|---|---|---|
| 14 | 2 mm grid; front speed off up to +50% (W22, `output/solver_verification.json`) | Fix | (a) Finish the running 1 mm pilot (27 pairs; interim: about 5 of 7 pairs per shard done at last check). (b) Pre-set rule: if 1 mm vs 2 mm Dice shift > 0.005 on the pilot, re-run the selected cells at 1 mm for all pairs on Kaggle GPU. | `output/pde_manifest/cache_h1/` | `output/pde_manifest/pilot_1mm.json` | Rule applied, result logged |
| 15 | Clamp adds mass, up to 42% (aniso r=10; W21) | Fix | Replace explicit stencil with a **positivity-preserving** scheme (non-negative directional splitting, or a monotone flux limiter). Log mass every run. Re-run solver tests V1–V9. If the r=10 arm still needs > 1% clamp mass, drop r=10 as an arm and say why. | Monotone schemes via non-negative directional splittings (Commun. Comput. Phys.) [S]; Lipnikov et al. monotone finite-volume schemes [S] | `src/solver_monotone.py`, `output/solver_verification_v2.json` | Clamp mass < 1e-6 on all tests |
| 16 | Grid edge, d = 0.003 in one fold (W3) | Fix | Extend d grid down one step (0.001) on development folds only, or use a bounded continuous optimiser. Report whether any fold still sits at an edge. | — | `selected_cells.json` (new run) | No fold at edge, or edge stated |
| 17 | D not identifiable (124 of 154 at lower bound; W14) | Finding + Fix | Stop reporting D and ρ alone. Report the identifiable pair: infiltration length λ = √(D/ρ) and front speed v = 2√(Dρ). Run synthetic recovery and profile likelihood. Use population priors (Bayesian). | Two-contour inference recovers λ well and v less well (Ezhov et al., arXiv 2111.13404) [S]; travelling-wave relations [S]; speed-only methods cannot separate D and ρ [S] | `output/identifiability.json` | Recovery plot + profile curves |
| 18 | Noise floor is an assumption (W24) | Fix | Measured floor from #12 (two segmenters + RHUH expert masks) and from literature repeatability (31–46%, AJNR) [S]. Recompute SNR at 1 mm. | as #4, #12 | `output/forecastability_measured.json` | Assumed and measured floors side by side |
| 19 | 2,000 bootstrap draws, plan says 10,000 (W25) | Fix | Re-run every primary and secondary contrast with 10,000 draws. Fixed seed. | — | all final JSONs | Every final JSON says `n_boot: 10000` |
| 20 | No patient-specific fit; "digital twin" not supported | Fix (test) | Patient-calibrated model for the 56 patients with 3+ scans (2p): fit on scans 1–2 with a population prior, forecast scan 3. Compare with population cell and persistence. Name stays "patient-calibrated model". Use "digital twin" only if it beats both, pre-declared. | Hierarchical Bayesian calibration (GPHMC, HAL inria-01324849) [S]; reduced-order + neural nets for patient-specific growth (arXiv 2412.05330) [S] | `output/patient_calibrated.json` | Rolling-origin result on 56 patients |

### D. Treatment and RL (simulation only)

| # | Blocker | Type | What we do | Literature / data | Output | Done when |
|---|---|---|---|---|---|---|
| 21 | Kill rates assumed, same for all 21 (W13) | Fix (try) / Accept | (a) Use the Rockne proliferation-invasion-radiotherapy model with linear-quadratic radiation term. In that work, radiation response correlated with ρ (r = 0.89, 9 patients) [S]. (b) Fit a per-patient response from MU on-treatment volume change where RT dates exist (24 patients lack dates, 2v). (c) Hierarchical prior. (d) If not estimable, keep "assumed" and move to supplement. | Rockne et al. 2010 Phys Med Biol 55:3271 [S] [VERIFY]; "Days Gained" metric, 33 patients, cut-offs 100 / 117 d (Neal et al. 2013 PLoS ONE) [S] [VERIFY] | `output/kill_rate_fit.json` | Either fitted per-patient values with CI, or a logged "not identifiable" |
| 22 | ρ > 0.02 cut-off chosen after data (W17, n = 10) | Fix | Replace with a continuous model (logistic: adaptive-earlier vs ρ) with a confidence band. Lock any threshold on training data. LUMIERE has only 4 and 3 patients above 0.02 (2l), so an external test is underpowered; say so. | — | `output/rho_continuous.json` | Band plotted; cut-off called "post hoc" everywhere |
| 23 | Equal-budget RL is simulation | Accept | Keep in supplement as "simulation". Add stress tests already planned (noise, misspecification). A one-line rule matches RL (C1), so RL is not a headline. | — | — | Wording fixed in paper |

### E. Molecular track (Track A)

| # | Blocker | Type | What we do | Literature / data | Output | Done when |
|---|---|---|---|---|---|---|
| 24 | 21 patients, SD up to 13 points (W9) | Fix | Replicate on GBmap (over 1.1 million cells, 240 patients) with patient-level splits. Leave-one-patient-out on the original 21 with patient bootstrap CI. | GBmap, Ruiz-Moreno et al. (PMC12526130) [S] | `output/trackA_gbmap.json` | Patient-level accuracy with CI on both atlases |
| 25 | "Healthy" = 3 non-GBM donors (W10) | Fix | Drop Healthy from all claims. Test Core vs Periphery within the 8 patients who have both. Add Ivy GAP anatomic regions (leading edge, infiltrating, cellular tumour) as a second region-labelled set. | Ivy GAP [VERIFY]; Ravi et al. 2022 Cancer Cell spatial programmes (PubMed 35700707, data on Dryad) [S] | `output/trackA_core_vs_periphery.json` | Within-patient test reported |
| 26 | Cells-per-patient shortcut 64.6% vs model 63.8% (2o) | Fix | Equal cells per patient (balanced subsample). Pseudo-bulk per patient per zone. Re-run all baselines. | — | `output/trackA_balanced.json` | Model vs shortcut on balanced data |
| 27 | Graph and cVAE leak (W12); MGMT underpowered (W18); zone files likely synthetic (W19) | Fix / Accept | Drop leaky models from claims (keep only `gat_pca_clean`-style). MGMT: report minimum detectable HR 1.74 (2j), context only. Zone files: never use; delete from any figure. | — | — | Removed from paper and board |

### F. Not done yet

| Item | Where in this plan | Phase |
|---|---|---|
| 1 mm pilot (running) | #14 | P0 |
| Perturbation test (`src/103_segmentation_perturbation.py`, written, not run) | #4, #18 | P1 |
| Negative controls | P1: shuffled patient pairs, time-reversed pairs, random-direction tensor, wrong-patient tensor | P1 |
| Compartment forecasts | P3: enhancing vs necrotic vs FLAIR vs cavity, scored separately | P3 |
| Freeze | #13 | P4 |
| LUMIERE external run | H-3, after freeze | P4 |
| Failure atlas | P4: 13 failure modes from masterplan; top 20 worst patients with images | P4 |
| `docs/isef_compliance.md` (author's fields) | P5. Author fills; Claude cannot fill personal fields | P5 |

---

## 4. Phases and dates

ISEF 2027 is May 8–14, 2027 (`paper/masterplan_implementation_plan.md`). Regional fair date: [author to fill].

| Phase | Dates | Work | Gate (must pass before next phase) |
|---|---|---|---|
| **P0 Setup + lock** | Oct 8–12 | Docker test; Kaggle notebook with GPU solver check; finish 1 mm pilot; 10,000-draw rerun of existing contrasts (#19); write and commit **Amendment 6** (H-2, H-3, rules below) | Amendment 6 committed with date, before PREDICT-GBM download |
| **P1 Numerics + noise** | Oct 13–26 | Monotone solver (#15); extended grid (#16); measured noise floor (#12, #18); perturbation test; negative controls; volume-matched shape test (#2); first-vs-later swap (#5) | Solver tests V1–V9 pass on new solver |
| **P2 H-2 recurrence** | Oct 27–Nov 23 | Download PREDICT-GBM; register tensor atlas; run iso and aniso models on TUM-GBM (development); freeze; run once on LUMIERE + RHUH part (test) | Freeze tag for H-2 set before test run |
| **P3 H-3 volume** | Nov 24–Dec 14 | Volume ladder on MU; conformal intervals; compartment forecasts; IPW sensitivity (#9) | Calibration plot exists |
| **P4 Patient fit + freeze + external** | Dec 15–Jan 11 | Patient-calibrated model (#20); identifiability (#17); kill-rate fit (#21); continuous ρ (#22); freeze `v1.0-development-frozen`; LUMIERE external run; failure atlas | Tag set before external run |
| **P5 Track A + paper** | Jan 12–Feb 15 | GBmap replication (#24–26); paper draft from ledger; board; ISEF forms | Every paper number cites a file |
| **P6 Review + practice** | Feb 16 → fair | ChatGPT review round 3; fix what checks out; judge drills | — |

Order if time runs short: P0 → P1 → P2 → P3. P4 and P5 can shrink to "exploratory" without hurting H-1, H-2, H-3.

---

## 5. Amendment 6 (draft content; commit before any H-2 or H-3 run)

H-2 primary: patient-level paired difference in **enhancing-recurrence coverage**, anisotropic model target minus standard 15 mm margin target, same volume, GBM test patients. Patient bootstrap 10,000 draws + Wilcoxon.
H-2 key secondary (Holm): anisotropic minus isotropic; anisotropic minus GliODIL and U-Net precomputed maps (same patient lists).
H-2 development set: TUM-GBM part (142). Test set: LUMIERE + RHUH parts. Tune only on development.
H-2 tensor: one pre-declared atlas; max 2 tensor variants (plain, free-water-corrected); Holm across them.
H-3 primary: mean |ln volume ratio| error, model minus persistence, patient level, MU development → LUMIERE frozen test.
H-3 secondary: 80% and 90% interval coverage; RANO-PD AUC on LUMIERE.
Abstain rule: SNR cut fixed on MU before LUMIERE.
Eligibility: pair-table rule (W27), 133 analysed (105 GBM).
Success words: "broadly positive" only if the 95% CI is above 0. Otherwise "not supported".

---

## 6. Decision rules (what if it fails)

| Result | What we say |
|---|---|
| H-2 positive (CI above 0 vs margin and vs isotropic) | "On a public 3-centre benchmark, tract direction improved recurrence coverage at equal volume." Main headline. |
| H-2 positive vs margin, null vs isotropic | "A growth model helps; tract direction does not." Matches MU. Still a clean answer to the benchmark's open question. |
| H-2 null | "Neither growth models nor tract direction beat the standard margin here." Pairs with H-1 into one story: limits of image-based GBM forecasting, measured. |
| H-3 positive | Volume forecast with calibrated uncertainty becomes the practical tool. |
| Everything null | Headline is H-1 with measured noise floor + 3 cohorts + pre-registration. Honest and strong for a fair: judges reward rigour. |

---

## 7. Words we never use unless a result earns them

"digital twin", "revolutionary", "first" (only after a literature check, cited), "validated" (only for external frozen tests), "predicts survival", "personalised therapy". Treatment results are always "simulation".

---

## 8. Literature found today (all [S]; check each before citing)

- PREDICT-GBM benchmark, 243–253 patients, 3 centres, recurrence coverage: arXiv 2509.13360 (https://arxiv.org/html/2509.13360v2); npj Digital Medicine (https://www.nature.com/articles/s41746-026-03194-0); code (https://github.com/BrainLesion/PredictGBM, Apache-2.0); data (https://huggingface.co/datasets/LZimmer/PREDICT-GBM, MIT, 25.5 GB).
- GliODIL, physics-informed radiotherapy planning, 152 patients: Nature Communications (https://www.nature.com/articles/s41467-025-60366-4).
- RHUH-GBM, 40 patients, pre-op / early post-op / recurrence, expert-corrected masks: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10551826/
- LUMIERE, 91 patients, 638 study dates, RANO ratings: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9755255/
- BraTS 2025 LUMIERE response task entry (AUC 0.81, macro F1 0.50; preprint): https://arxiv.org/abs/2509.06511
- Deep learning for reaction-diffusion parameters; λ recovered better than v: https://arxiv.org/pdf/2111.13404
- TaDiff, treatment-aware diffusion model for glioma growth: https://pubmed.ncbi.nlm.nih.gov/40031286/
- Multi-task diffusion for glioma progression (preprint): https://arxiv.org/pdf/2509.10824
- Rockne et al. 2010, radiotherapy response model: https://mayoclinic.elsevierpure.com/en/publications/predicting-the-efficacy-of-radiotherapy-in-individual-glioblastom
- Jackson et al. 2015 review of proliferation-invasion model uses: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4445762/
- Unkelbach et al., RT planning with a growth model, 10 patients: https://arxiv.org/pdf/1311.5902
- BraTumIA test-retest repeatability: https://www.ajnr.org/content/42/6/1080
- Automated volumetry ICC 0.97 vs manual 1D/2D: https://arxiv.org/pdf/2209.01402
- Free-water-corrected FA and recurrence: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12696142/ ; FA + free-water recurrence prediction: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7140058/
- DTI for GBM progression, systematic review: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11626419/
- Conformal volume under shift (MICCAI 2024): https://arxiv.org/pdf/2407.19938 ; conformal bands for irregular visits: https://arxiv.org/pdf/2511.13911 ; sparse lesion forecasting (preprint): https://arxiv.org/pdf/2609.21197
- Informative dropout weighting: https://arxiv.org/pdf/2510.19154
- Pre-op growth rate and survival, mixed evidence: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11046984/ ; https://ar.iiarjournals.org/content/44/11/5043
- Monotone anisotropic diffusion schemes: https://core-cms.cambridgecore.org/core/journals/communications-in-computational-physics/article/monotone-finite-difference-schemes-for-anisotropic-diffusion-problems-via-nonnegative-directional-splittings/A380E9FC615239815892D6D360EEA310
- GBmap: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12526130/ ; Ravi et al. 2022: https://pubmed.ncbi.nlm.nih.gov/35700707/
- UPenn-GBM (630 patients, has DTI, mostly pre-op): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9338035/

Not found in today's search: a published GBM study that compares a PDE forecast with a persistence baseline on next-scan shape. If that holds after a proper search, H-1 itself is a gap filled. [VERIFY with a PubMed search before claiming.]

## 9. Limit removal (added 2026-10-08)

Rule from the author: find a fix for every stated limit. Rule from CLAUDE.md: a limit that is not fixed stays in the paper.
Each row: limit, fix, source, status. "Needs OK" = needs the author's permission (download or upload).

| # | Limit | Fix | Source | Status |
|---|---|---|---|---|
| L1 | Monotone solver not exact where a tensor is not a non-negative sum of 37 offsets (99.66% exact at r = 10) | Selling's decomposition: 6 non-negative weights, exact for every SPD tensor; long offsets blocked if they cross non-domain voxels | Selling 1874; Conway & Sloane 1992; Fehrenbach & Mirebeau 2014 (arXiv 1301.3925) | Done in code (`split="selling"`, 8 tests pass); script 97 rerun running |
| L2 | V7: 2 mm grid is coarse against the front width (W22) | (a) 1 mm production run on Kaggle GPU; (b) discrete-speed correction: calibrate D_eff(h) so the discrete front speed equals 2 sqrt(D rho), applied before the fit; (c) anisotropic eikonal forecast arm, which does not need to resolve the front width | (c) Konukoglu et al. 2010, IEEE TMI 29:77 | (a) Needs OK (private Kaggle dataset of derived 2 mm masks); (b) can run locally; (c) Amendment 7 before any test-set use |
| L3 | Production aniso arms not rerun with the new solver | Rerun script 100 aniso arms with `TensorFKMonotone` (Selling is about 8x slower at r = 10 on CPU, so GPU) | — | Open; GPU preferred |
| L4 | Script 103 noise size assumed (1 voxel) | Script 106: measured boundary disagreement (Dice, ASSD, HD95) between two segmenters on all LUMIERE scans, native and 2 mm | LUMIERE native label maps | Running (declared in 4a8091e before the run) |
| L5 | Two automated tools are not a scan-rescan pair | RIDER Neuro MRI was the candidate; since April 2025 it is controlled access via dbGaP (needs a PI, eRA Commons, institutional signing official) | TCIA Legal Framework page | Accepted (author decision 2026-10-08): paper states "No true scan-rescan data was available to us"; noise is bounded by scripts 103, 104, 106 |
| L6 | LUMIERE enhancing label is not MU core (HD-GLIO-AUTO has no necrosis label) | Run one segmenter with the same labels on MU and LUMIERE | BraTS-trained model | Partial: needs GPU and a license check before any LUMIERE upload (non-commercial licence) |
| L7 | Noise floor per volume tertile, not continuous | Smooth fit of disagreement against ln volume; sensitivity next to the registered A6.5 floor | — | Can run locally (sensitivity, not primary) |
| L8 | Scans where one tool finds nothing are excluded | Include with a +1 voxel offset, as script 99 does; sensitivity | — | Can run locally |
| L9 | Script 105: PDE density only inside a box around the input core | Larger box / whole-brain domain (GPU) | — | Open |
| L10 | One population cell per fold; no patient-specific fit | Fit each patient's cell on the previous pair, forecast the next | — | Open |
| L11 | Retrospective only; automated masks are not ground truth | Cannot be removed by analysis | — | Accepted (W20); stated |


## 10. Gap analysis and fresh literature (merged from GAP_ANALYSIS.md, 2026-10-09)

End goal: (H-1) measure how forecastable glioblastoma growth is; (H-2) test whether a growth model can place the
radiotherapy margin better than the guideline 15 mm margin, at the same treated volume; (H-3) forecast next-scan volume
with calibrated intervals and an abstain rule. Numbers for our results: `vault/PAPER_FINDINGS_LEDGER.md`.
Literature below was found by web search on 2026-10-09; each entry is [S] (search summary) until the full text is read.

### 10.1 Where the field is (fresh literature)

| Finding | Source | What it means for us |
|---|---|---|
| Guideline margin is 15 mm around enhancing tumour / cavity (was 20 mm in 2016); no consensus on a T2/FLAIR margin (0-15 mm) | Niyazi et al. 2023, ESTRO-EANO, Radiother Oncol 184:109663 | M0 in H-2 = this plan; confirm PREDICT-GBM `create_standard_plan` uses the same rules |
| PREDICT-GBM: 243 patients, multi-centre; at equal treated volume, U-Net 79.37 +/- 2.08% and GliODIL 78.91 +/- 2.08% enhancing-recurrence coverage; only U-Net and GliODIL significantly beat the standard plan pooled; LMI and nnU-Net were below it; distant recurrences missed by all; gains "modest" | Zimmer et al., arXiv 2509.13360 (v2 Mar 2026); npj Digit Med 2026 | The bar to beat is about 79%. A plain population Fisher-Kolmogorov (like LMI) can fall BELOW the standard plan. Our M1/M2 need per-patient fitting to compete |
| GliODIL: physics-informed discrete loss, infers full tumour-cell map from MRI (+ FET-PET); best recurrence prediction among tested models, 152 test patients | Balcerak et al. 2025, Nat Commun, doi 10.1038/s41467-025-60366-4 | The strongest biophysical competitor; it fits per patient. Its predictions are in the PREDICT-GBM release, so we can compare without running it |
| Prospective pilot: ML-guided personalised dose escalation in new GBM; survival comparison preliminary | Nat Commun 2026, s41467-026-72545-y | Field is moving to clinic; context for the poster, not a method to copy |
| Atlas fibre orientation warped to the patient gives a significant but small Dice gain, mostly in butterfly gliomas crossing the corpus callosum | arXiv 2507.17707 (2025) | Matches our null overall (ledger 2ak). Suggests a pre-registered subgroup: tumours touching the corpus callosum |
| No head-to-head study of patient DTI vs atlas tensor for recurrence | search 2026-10-09 | Open question we could answer with UPENN-GBM DTI (317 patients) if recurrence or follow-up masks exist |
| No shared benchmark for next-scan glioma prediction; treatment-aware diffusion model predicts future masks with uncertainty | IEEE TMI 2025 (PubMed 40031286); arXiv 2509.10824 | H-3 has no standard competitor; persistence and the geometric rule stay the baselines |
| Conformal intervals for longitudinal lesion size; time-varying nonconformity score for irregular follow-up; exchangeability breaks under site shift (weighted conformal) | arXiv 2609.21197 (2026); Tassopoulou et al. (OpenReview 2025); arXiv 2407.19938 | H-3 method: split conformal on MU, time-varying score (our intervals are irregular), report coverage on LUMIERE as a shifted site |

### 10.2 What we are missing

#### H-2 (the possible strong positive)
1. **PREDICT-GBM data not downloaded** (Hugging Face LZimmer/PREDICT-GBM, about 25.5 GB). Needs the author's OK. Amendment 6 is already committed, so the download is allowed by our own rules.
2. **Per-patient fitting.** Our models use one population cell per fold. The PREDICT-GBM result says population models (LMI) can lose to the 15 mm margin. Fix: fit lambda = sqrt(D/rho) per patient from the pre-op scan (enhancing edge vs FLAIR edge; Konukoglu 2010 eikonal idea; Lipkova 2019 Bayesian) using the exact solver. Must be added to Amendment 6 as a model before any test-set run.
3. **Adapter to PREDICT-GBM inputs** (their tissue maps, SRI24/MNI space, our atlas tensor registered to that space) and use of `topk_plan` for equal volume.
4. **Compare against published predictions** (U-Net, GliODIL, LMI files in the release): report our model next to them on the same patients, not only M2 vs M0.
5. **Distant recurrence**: nobody captures it. Report the share of distant recurrence separately; do not hide it.

#### H-3
6. **Conformal intervals and abstain rule not built.** Split conformal on MU development folds, time-varying nonconformity score, then frozen test on LUMIERE; report coverage at 80/90%.
7. **LUMIERE volumes are noisy** (tools disagree about follow-up change by more than the change, ledger 2ab). State that the LUMIERE test bounds what any forecast can show.

#### Direction question (follow-up to ledger 2ak)
8. **Patient DTI test.** UPENN-GBM has DTI for 317 patients (CC BY 4.0). Check first whether it has follow-up or recurrence scans; if not, patient DTI can only be tested where both DTI and a later mask exist.
9. **Corpus-callosum subgroup**, pre-registered before looking.

#### H-1
10. **True scan-rescan noise**: not available openly (RIDER Neuro MRI and QIN-GBM are dbGaP controlled access). Measured substitutes: tool disagreement (2ab, 2ag), expert correction (2ai). Simulated rescans (scanner perturbation, then re-segmentation) are still open.

#### Paper and poster
11. Remove every "direction helps" claim (ledger 2ak shows it was a solver artefact).
12. Write the solver-artefact finding as its own result, with the Selling method and the GPU rerun.

### 10.3 Order of work (recommended)
1. Fix old claims (item 11). Small, urgent.
2. Amendment 7: per-patient fitting model and the corpus-callosum subgroup, committed BEFORE the PREDICT-GBM download.
3. PREDICT-GBM download (needs OK), adapter, development-set runs on TUM only.
4. Freeze (`configs/h2_frozen.yaml` + tag `h2-frozen`), then test on LUMIERE + RHUH.
5. H-3 conformal (CPU, can run in parallel).
6. UPENN-GBM DTI check (item 8).

### 10.4 Sources
- Niyazi et al. 2023 ESTRO-EANO: https://pubmed.ncbi.nlm.nih.gov/37059335/
- PREDICT-GBM: https://arxiv.org/abs/2509.13360 ; https://www.nature.com/articles/s41746-026-03194-0
- GliODIL: https://www.nature.com/articles/s41467-025-60366-4
- ML-guided dose escalation pilot: https://www.nature.com/articles/s41467-026-72545-y
- Butterfly glioma fibre-tract model: https://arxiv.org/pdf/2507.17707
- Treatment-aware diffusion model: https://pubmed.ncbi.nlm.nih.gov/40031286/
- Multi-task diffusion glioma progression: https://arxiv.org/html/2509.10824v1
- Conformal lesion-size forecasting: https://arxiv.org/pdf/2609.21197
- Robust conformal volume estimation: https://arxiv.org/pdf/2407.19938
- Konukoglu et al. 2010 IEEE TMI: https://nmr.mgh.harvard.edu/node/4758


## 11. The final step: a 3D tumour forecast you can see (added 2026-10-09)

What it is: an interactive 3D page (rotate, zoom, switch layers) for one patient, built from our real outputs:
1. Brain outline (grey, see-through).
2. Tumour today (input core mask).
3. Model forecast for the next scan (PDE density as a colour cloud, plus the visible-tumour threshold surface).
4. What really happened at the next scan (true mask), so the viewer sees where the model was right and wrong.
5. For H-2 patients: the standard 15 mm radiation volume vs the model-shaped volume of the same size, and the real
   recurrence. This is the picture that shows the whole project in one view.

How: marching cubes (scikit-image) turns each 3D mask into a surface; the page draws the surfaces in the browser
(three.js or Plotly mesh3d), one checkbox per layer, a slider over forecast days. Script number: next free (110).
Data rule: the page uses masks only (no MRI intensities); it stays private unless the author decides to share it.
Honesty rule: the page shows the Dice for that patient and says "forecast", "atlas tensor", "retrospective".
Done when: the page opens from one file, works for any patient id, and the shown Dice equals the ledger/cache value.
Order: build it now on MU (we have forecasts and true next scans); add the radiation-plan layers after H-2 runs.
