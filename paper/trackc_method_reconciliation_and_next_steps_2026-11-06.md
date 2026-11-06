# Integrated GBM method review and next steps (2026-11-06)

**Status:** corrective methods audit in progress. The research title is unchanged: **“Integrating Tumor Biology and Longitudinal MRI for Glioblastoma Progression Modeling and In-Silico Therapy Optimization.”** No new policy-outcome comparison, patient analysis, or sealed evaluation was run for this note. Nothing here can establish that the project is the best, most unique, or guaranteed to produce a positive result.

## What prior work means for the claim

The broad combination in the title is already represented in prior work. The 2025 M4RL paper combines a multiscale glioblastoma tumor/macrophage model, a Fokker–Planck/PINN surrogate, and reinforcement learning to schedule CSF1R and IGF1R inhibitors ([Liu et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC13155510/)). A 2026 preprint combines reaction–diffusion growth, a learned residual, online digital-twin updating, and MPC treatment scheduling on 387 synthetic trajectories ([Liu et al., 2026](https://arxiv.org/abs/2607.13877)). These papers do not settle whether the project's narrower decision comparison is new, but they rule out presenting “biology + MRI + tumor model + adaptive therapy” alone as the novelty.

The defensible candidate question remains narrower: **under explicitly simulated links between paired baseline biology-like and MRI-like signals and latent resistance, when does dependence-aware policy updating help or harm compared with a matched conditional-independence policy?** A novelty claim for that intersection remains unverified and would require a focused review before publication.

## Corrections made to align plan, code, and evidence

1. **Utility estimand:** the code uses patient-level worst-endpoint utility before averaging within a world:

   `U_i = min(min(T_baseline_i, 730), min(T_RANO_i, 730)) - 5 * courses_i_through_day_730`.

   The previous plan wording could instead mean “take the smaller of the two cohort-average RMSTs.” Those are different estimands. The full and short plans now state the implemented patient-level definition. A regression example with endpoint pairs `(100, 700)` and `(700, 100)` gives a mean patient-level minimum of 100 days, whereas the minimum of endpoint means is 400 days.

2. **Separate synthetic cohorts:** the policy harness now requires unique, canonical patient IDs for selection and evaluation, rejects overlap before policy fitting, and fingerprints both ordered ID lists. It also records separately declared scenario-generation seeds and notes that they are caller-declared; hashes identify the actual scenario arrays consumed.

3. **Censoring audit:** checked arm outputs now retain both endpoint-cap flags. These flags are needed to distinguish a progression event from follow-up ending at the simulator horizon.

4. **Treatment paths:** policies share the same virtual-patient inputs, precomputed observation-noise stream, scan schedule, and treatment budget. Each policy has its own action-dependent disease trajectory. The two endpoint times for a given policy are scored on that policy's trajectory and course count. Descriptions that implied one trajectory shared across different policies have been corrected.

## Existing signal-prediction result that constrains the next method

The frozen leave-one-out (LOO) signal-prediction gate did not pass at cohort sizes 18 or 28. The LOO model average beat the uniform-grid mixture, but it lost to the matched conditional-independence comparator: mean log-score differences were −0.019283 (one-sided 98.75% lower bound −0.029323) at n=18 and −0.008081 (lower bound −0.014564) at n=28. The four-bound advancement rule therefore failed. The comparator was given the generating `q`, `a`, and `b`, while the LOO approach had to infer a shared model; this makes the comparison intentionally demanding, but does not erase the failed result. The LOO method will not advance as the primary policy method from that seed block. See the [frozen predictive-check protocol](trackc_loo_predictive_check_protocol_2026-11-06.md), [method note](trackc_loo_model_averaging_note_2026-11-06.md), and [saved diagnostic](../output/trackc_loo_predictive_check/diagnostic.json).

## Next method-development path

The next candidate should address the documented weak finite-sample identification of dependence rather than tune the failed LOO method on spent seeds. A separate, explicitly **proposed** protocol could:

1. Estimate a set or distribution of plausible signal-dependence parameters using development/calibration worlds only, retaining weakly identified values instead of selecting one point estimate.
2. Freeze a robust policy that evaluates actions across that uncertainty set, with no fitting or peer-weight adaptation on evaluation worlds.
3. Compare it on fresh, disjoint virtual worlds with a conditional-independence policy that receives the same calibration information, plus single-signal, no-signal, fixed-schedule, and equal-course controls.
4. Check signal prediction and policy value as separate gates. If signal prediction is inadequate, label policy findings as a simulator stress test rather than evidence that the estimated signal model works.
5. Freeze the world-level estimand, sample-size/precision target, bootstrap unit, simultaneous interval method, numerical tolerances, and all evaluation anchors before generating evaluation outcomes. Report all anchors, including harm; do not change the frozen rule to obtain a positive result.

This is a design direction, not an adopted method or a claim of superiority. A repository search found no research-scale driver calling `run_shared_policy_world`; current execution evidence is a small contract test, not a multi-world result. The harness now has explicit six-course, twelve-course front-loaded, and twelve-course 56-day pacing references. The 56-day arm is not a dose-dense regimen; the simulator's 28-day cycle cannot represent intervals shorter than 28 days. Before any new outcome comparison, an executable driver still needs a full audit for an equal-course fixed-schedule frontier, held-out world aggregation, a precision target, an interval method, and exact evaluation-anchor enforcement. The synthetic contract tests below establish plumbing only; they do not show a treatment-policy advantage.

## Current project boundary

Track A and the Duke patient analysis remain separate from Track C. There is no patient-level match between their cohorts, and the molecular-to-resistance relationship remains an explicit simulator assumption. The project will not describe simulated drug suppression as a real treatment recommendation or as evidence of clinical benefit. The title remains unchanged; any later plan amendment must be explicit and dated.
