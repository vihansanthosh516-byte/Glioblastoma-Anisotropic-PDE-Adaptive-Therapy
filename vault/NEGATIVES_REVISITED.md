# Six Negatives Revisited

**Date:** 2026-10-01
**Rule followed:** no committed result was changed except script 42's output, which was regenerated on request. Every new analysis is a new numbered script with its own JSON, and each one is committed separately. Primary comparisons were fixed before each run. Where one was changed after seeing data, that is stated.

| # | Negative | Correct test | Effort | P(positive) | Outcome | Recommendation |
|---|---|---|---|---|---|---|
| B1 | Forecast: no model beats no-change | Sub-region (FLAIR vs enhancing) scoring | 6-8 h (4.5 h compute) | Low | Not run; data analysis below | Document as data limitation |
| B2 | Ablation < 0.01% | Spatial endpoints (script 77); horizon x D x grid sweep (script 78) | 1 h + 21 h compute | Low -> High | **Parameter artifact, not structural.** Script 60 D is ~10x below the literature; at D = 0.1 the DTI effect is 8% Dice at 90 d, same at 1 mm | Report script 78; re-check script 60 at D = 0.1-0.13 |
| B3 | Script 42 blocked; D_f aniso vs iso | Matched-isotropic paired test (script 74) | 2 h | Medium | D_f: significant, wrong direction, not a valid metric. **Elongation: strongly positive** | Report elongation; drop D_f |
| C1 | Unconditioned PPO 25% | Probe-then-commit, kill rates inferred (script 75) | 2 h | Medium | **Positive noise-free (100%), fails at realistic noise** | Report both; limitation |
| C2 | Drug-budget artifact 100% -> 10% | Paced heuristic, day-90 and TTP (script 75) | 1 h | High | **Positive on TTP in 4/4 sets; negative on day-90** | Report: pacing extends TTP |
| C3 | No adaptive structure | Sensitive + resistant LV model + same-window controls (script 76) | 4 h (+12 h compute) | Medium | **Stays negative for real patients; resistance-driven gain only in 2/384 synthetic rows** | Document as limitation; do not retrain RL |

---

## B1 - Forecast (152 MU-Glioma-Post patients): data limitation

No new run. Evidence from the committed `output/forecast_validation/forecast_results.json`:

- **80/152 (53%) tumours shrank or stayed the same between scans.** No growth model can forecast shrinkage. The loss to no-change is concentrated there: aniso - no_change = -0.068 Dice (p = 1.4e-10) in shrank/same, against -0.023 (CI crosses 0) in grew, where aniso beats no-change in 61% of patients.
- **Direction adds nothing even when volume is known.** Volume-matched, growing only (n = 72): every FK arm beats uniform dilation (+0.0068 Dice, p = 0.017), so PDE-shaped growth helps a little. Aniso vs iso_same = -0.00003 (p = 0.65), so the DTI orientation itself contributes nothing.
- 117/152 had radiation and 108/152 had TMZ schedules. The treatment kill term (commit 9a9526a) made the forecast worse (b1fb529).

Options:

| Option | Literature | Data supports? | Verdict |
|---|---|---|---|
| Log-Euclidean tensor atlas | Arsigny 2006; standard DTI-atlas practice | Atlas FA is already lower than per-patient FA (dispersion QC), but aniso-vs-iso is 0.0002 Dice even with sharpening r on the grid | Low; would refine an effect of ~0 |
| Patient-specific DTI (UCSF-PDGM) | Jbabdi 2005; Painter & Hillen 2013 | One scan per patient, so no forecast target | Not feasible |
| Treatment kill term | - | Already done; regressed | Done, negative |
| Sub-regions (FLAIR vs enhancing) | Swanson 2008 thresholds (T2 ~16%, T1Gd ~80% of max density); anisotropic invasion appears on T2/FLAIR | Masks carry labels (`*_tumorMask.nii.gz`, currently `> 0`) | Best remaining option, but low probability. Shrinkage dominates whole-tumour change and would dominate sub-regions too |

## B2 - Ablation: parameter artifact, not structural (scripts 77 and 78)

Full write-up: [Script-78-Horizon-Crossover.md](Script-78-Horizon-Crossover.md).

- Script 60 scores total volume ∑u. Diffusion conserves mass, so volume is insensitive to D by construction. Script 77 confirmed that spatial endpoints (mask Dice 0.994, extent) also show no effect at script 60's settings.
- Script 77 attributed that to the length scale: over 90 d the diffusion length is 0.46-1.18 mm against a 2 mm voxel. That is correct **for D_white = 0.013 mm²/d, which is ~10x below the Swanson-type 0.13** (D_gray 0.0013 vs 0.013). Script 60's values match the cm²/d numbers read as mm²/d.
- Script 78 (D x horizon x grid sweep, 16 cells, endpoints fixed in the header):
  - Dice(full vs matched-isotropic, 2 mm, 90 d): 0.997 (D 0.01), 0.985 (0.03), **0.923 (0.1), 0.838 (0.3)**.
  - Crossover (Dice < 0.95) at 90 d for D ≥ 0.1; none for D ≤ 0.03.
  - 1 mm grid gives the same effect (0.943 at D 0.1, 0.848 at D 0.3), so it is physics, not resolution.
  - The effect is largest early and shrinks at 365-900 d (the tumour fills the domain).
- Caveats: modest effect (Dice 0.92-0.95 at D = 0.1), no noise comparison, rho 0.01 cells censored at 90 d (mask under threshold), elongation not monotone. Details in the script-78 note.
- Verdict: the earlier "structural" label was wrong. Anything concluded from script 60 about DTI (including Sobol S1(rho) = 0.998) used a diffusivity 10x too low.

## B3 - Script 42 + anisotropic vs isotropic (commits b2a395a, 109e0d4)

- **Script 42 fixes:**
  - Resolved the 6 conflict blocks and fixed `base=2`.
  - Phase 3 had simulated **0 patients**: the location CSV uses `PAT_xxxx` IDs and the NPZ uses `PatientID_xxxx`. It now simulates all 103.
- **Script 42 results:**
  - D_f median 0.909, with 0/103 > 1.2.
  - **D_f is not in 1.0-2.0 on either side**: aniso 101/103 < 1, iso 100/103 < 1. The tumours are ~120-140 px, so box counting saturates and measures size.
- **Script 74, paired test, n = 103.** Only the tensor differs between arms: iso = tr(D)/2·I.
  - D_f: aniso 0.925 vs iso 0.943, diff **-0.018, p = 2.9e-15**. Significant, but the **opposite** of the hypothesis, and driven by area (aniso 123 vs iso 137 px). D_f is not a valid anisotropy metric here.
  - **Elongation: 1.21 vs 1.03, diff +0.178, d_z = 2.9, p = 4.1e-19, 98/103.** At 600 d: 1.56 vs 1.13.
  - **Positive finding: anisotropy produces a large, consistent shape effect, detectable by elongation, not by D_f.**
  - Script 45's isotropic cache (D_f 1.009, 8 legacy patients, different solver settings) is not a fair comparator.

## C1 - PPO without kill rates (script 75, commit 25c027a)

Blind policies fail because the day-90 optimum is a late block, so nothing informative is observed before the action must be chosen. The fix is to pay for information with a probe:
1. Give off/chemo/RT for L days each.
2. Infer k = g_off - g_drug from the volume response.
3. Commit using the inferred kills.

| cohort64 (primary: L=1, noise-free, rule) | day-90 win vs Stupp | TTP vs Stupp |
|---|---|---|
| PPO66 (blind) | 25% | +34.7 d |
| script-67 blind nets | 40-45% | 0 d |
| **probe + rule, sigma 0** | **100%** (log ratio -0.476 [-0.663, -0.322], p = 1.9e-6) | +2.6 d |
| probe + rule, sigma 0.02 / 0.05 / 0.10 (3-seed mean) | 72% / 63% / 62% | +2.6 d |
| probe + TTP net, sigma 0 / 0.10 | 70% / 60% | **+58.9 / +50.4 d** |

- Under noise, lhs60 and real_test drop to 11-46%, worse than the blind nets there (100%). Kill rates are identifiable from daily noise-free volumes, not from clinical volumetry (~5-10% CV).
- **Verdict:** PPO can be fixed without being given kill rates only if response is measured precisely. Report this as an identifiability result, not a deployable policy.
- **Exception:** the TTP-trained net with a probe keeps +46-58 d on cohort64 at sigma 0.10. On real_test it loses (-20 to -35 d), so it is set-dependent.

## C2 - Paced heuristic at equal budget (script 75, commit 25c027a)

Cumulative AUC is capped at budget × day/90.

| set | day-90 win | TTP vs Stupp (95% CI) | longer/equal/shorter |
|---|---|---|---|
| cohort64 (20) | 0% | **+48.5 d [34.9, 61.7]** | 20/0/0 |
| lhs60 (30) | 0% | **+52.5 d [39.0, 66.1]** | 30/0/0 |
| real_test (21) | 0% | **+2.1 d [1.2, 3.3]** | 15/6/0 |
| synth_test (40) | 15% | **+15.7 d [7.8, 25.7]** | 36/4/0 |

- **Positive on the clinical endpoint**: TTP is longer in 101/111 patients and shorter in none, at equal budget, all p < 0.001.
- It loses on day-90 volume, which is the endpoint script 68 already showed rewards back-loading.
- The paced arm amounts to an evenly spread RT/chemo schedule. The volume threshold does nothing, so the finding is "pacing", not "the heuristic".

## C3 - Resistance dynamics (script 76)

**Model.** Sensitive and resistant cells share the carrying capacity K (Lotka-Volterra, as in Strobl 2021 and Gallagher 2024). Each set keeps its own rho and kill table, and the Stupp budget applies. With f_r0 = 0 the model reproduces the original exactly.

**Grid.** f_r0 {0.001, 0.01, 0.1} × fitness cost {0, 0.25} × residual kill on R {0, 0.2} × occupancy {script-59 seed, near capacity} × window {90, 365}. Plus 4 single-population controls. 4 sets, 111 patients, RANO TTP to day 365.

**Confound caught.** On the first run, AT80 beat Stupp by about 80 d with a 365-day window, even with almost no resistant cells. The cause was spreading the budget over a year, not adaptivity: paced 347.8 vs AT80 346.5 d on real_test. Adaptive arms are therefore compared against:
- a **same-window paced schedule**;
- the **per-patient best non-adaptive arm**;
- the **f_r0 = 0 control**, giving a difference-in-differences.

The primary window was set to 365 after that run. That change is disclosed in the script and commit.

| Primary cell (f_r0 0.01, cost 0.25, 365 d) | AT80 - paced (95% CI) | gain in f_r0=0 control | resistance-attributable |
|---|---|---|---|
| cohort64 | +90.5 [78.0, 104.0] | +93.3 | **-2.9** |
| lhs60 | +34.5 [-21.9, 84.0] | - | - |
| real_test | -19.7 [-49.7, 1.8] | - | - |
| synth_test | +11.8 [-7.6, 31.4] | - | - |

- **Adaptive beats the best non-adaptive arm in 19 of 384 cell × set × arm rows.**
  - 17 of those are cohort64, where the gain is identical without resistance. That is volume feedback on that set's kill rates, not clonal competition.
  - Only **2 rows** (lhs60, f_r0 = 0.1, near capacity or with a cost) show a resistance-attributable win: +22.3 and +15.9 d vs best non-adaptive. This is exactly the corner Strobl 2021 predicts.
  - **real_test: 0 rows.** AT50 loses everywhere.
- **Why:** the median resistant fraction at day 365 is about 2%. At these growth rates resistance barely expands within a year, so progression is sensitive-cell regrowth after treatment stops. There is little competition for adaptive timing to manage.
- **Verdict:** adding resistance does not turn this negative positive for the real-patient sets. Retraining the RL (scripts 66-68) on this model is not justified. Report the 2-row corner as consistent with theory.

---

Sources: Leder 2014 Cell 156:603; Gallagher 2024 Cancer Res (10.1158/0008-5472.CAN-23-2040); Strobl 2021 Cancer Res (turnover / cost of resistance); Zhang 2017 Nat Commun 8:1816; Swanson 2008 (T2 16% / T1Gd 80%); Painter & Hillen 2013 J Theor Biol; Jbabdi 2005 MRM; Arsigny 2006 MRM.
