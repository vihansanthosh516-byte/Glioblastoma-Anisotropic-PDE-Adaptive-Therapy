# Research Plan Draft — Track C: Robust Belief-Based Dose Timing in a Simulated GBM Model

Status: DRAFT (2026-10-27, Claude). Track C only. This file does not replace `paper/research_plan_complete_project_form.md`; the Track C section there is the earlier plan and is superseded only if the author accepts this draft. Every number below cites an output file. All results are computer simulation with assumed parameters. Nothing here is a patient, drug, or clinical result.

## Rationale

Adaptive therapy times drug courses to tumor burden instead of giving the maximum dose up front (Zhang et al., 2017; West et al., 2023). Model-based schedules are often tested only inside the model that produced them. A schedule can then look good because it assumes the right drug strength, growth rate and noise. This project asks a narrower question: does a dose-timing rule that acts on a *probability of near-term progression* keep its benefit when the simulator's truth departs from the controller's model?

Earlier exploratory and sealed work in this repository (disclosed in full below) found one rule, "wide-prior containment", that gained time to progression in the nominal model and showed no significant harm in 17 pre-set wrong-model tests. That finding has known weak points: the tests were chosen by the author, one test family (weak drug) was seen during earlier tuning, "no harm" was defined only as "no confidence interval entirely below zero", and two suites had wide intervals. This plan is a new, sealed confirmation designed to address those weak points. It is not assumed to succeed.

Novelty: not claimed. Risk- or probability-triggered dosing, particle filtering for tumor state, and robustness checks of adaptive therapy all have prior work. A focused prior-art search (terms: chance-constrained adaptive therapy, belief-space control cancer, particle filter adaptive therapy, robust adaptive therapy misspecification) is required before any novelty statement. Not done yet.

## What already exists (disclosed; all simulation, all exploratory except where marked sealed)

| Script | What | Result | Weakest point | Output |
|---|---|---|---|---|
| 150 | Pilot 2, original CVaR-tuned filter hybrid | Lost to front-loaded dosing on CVaR_0.10: 152.4 vs 196.1 days, diff −43.7 (CI −62.6 to −27.2) | Pilot seeds | `output/trackc_pilot2/pilot_results.json` |
| 153 | Known-parameter oracle headroom | Oracle beats best rule by mean 83.5 days; CVaR_0.10 headroom only +3.4 (CI 1.3 to 6.9) | Oracle is best-found, not proven optimal | `output/trackc_headroom/headroom_results.json` |
| 154 | Containment rule (narrow prior) | +33.0 days mean TTP vs raw on/off 0.6/0.5 (CI 26.8 to 39.1) | Endpoint switched to mean TTP **post hoc** | `output/trackc_containment/containment_results.json` |
| 155 | Wrong-model suites, narrow | Weak drug −76.8 (CI −86.4 to −66.6) | Revealed the failure mode | `output/trackc_wrong_models/wrong_model_results.json` |
| 156 | Wide prior on drug strength m ∈ [3, 20] | Weak drug +0.4 (CI −1.7 to 2.7); nominal +16.0 | Weak-drug range used in tuning | `output/trackc_guarded/guarded_results.json` |
| 166 | Filter calibration, **sealed** | Predictive 90% coverage 0.89; resistant-fraction coverage 0.73–0.77 (overconfident); high resistance not detectable from readings | Calibration of latent resistance fails | `output/trackc_filter_calibration/sealed_results.json` |
| 167 | Equal-budget tuning, **sealed** nominal test | Narrow +24.7 (CI 19.6 to 29.8), rule met; wide +21.1 (CI 16.0 to 25.9); front-loaded best CVaR (182.2 vs 177.7) but mean −13.1 | Hypothesis came from post hoc work; nominal model only | `output/trackc_sealed_containment/sealed_results.json` |
| 168 | Wrong-model grid, 17 suites, **sealed** | Wide: no suite with CI entirely < 0; CI entirely > 0 in 14/17; worst weak_drug_sev +1.2 (CI −3.1 to 5.4). Narrow: weak_drug_sev −81.7 (CI −98.1 to −65.5) | Author-chosen suites; n 150 too small for noisy suites (noise_x3 CI −17.9 to 29.0; combo_bias_noise CI −9.7 to 22.2) | `output/trackc_sealed_wrong_models/sealed_results.json` |

Important negative facts kept in view:
- Under a stricter non-inferiority margin of −7 days (used below), script 168 wide containment would **fail** in 2 of 17 suites (noise_x3 lower bound −17.9; combo_bias_noise lower bound −9.7). The 168 "no harm" result is absence of detected harm, not shown non-inferiority.
- Front-loaded dosing beats containment on the lower tail (CVaR_0.10) in the nominal model and under fast growth (168: fast_growth_mod +47.2 vs +33.2). Containment is a mean-TTP rule, not a tail-risk winner.
- The filter is overconfident about the hidden resistant fraction (166). No claim is made that the rule estimates resistance well.
- Per-scenario filter outcomes are Monte Carlo-sensitive: 23–25% of scenarios move > 1 day between step sizes (`output/trackc_numerics_check/numerics_results.json`). Only averages are interpreted.

## Research Question and Hypotheses

**Question:** In a two-population (sensitive/resistant) GBM growth simulator with a 12-course budget, does wide-prior containment remain non-inferior to the strongest equally tuned non-filter rule when the truth departs from the controller's model in ways **not used in any earlier tuning or testing**, and does it keep a mean gain when model departures are sampled at random instead of chosen one at a time?

**H1 (primary, robustness):** In every suite of a new sealed misspecification grid, the 95% lower confidence bound of mean TTP (wide containment minus comparator) is above −7 days.

**H2 (co-primary, sampled uncertainty):** In a suite where each virtual patient's departures from the model are drawn at random from a pre-set mixture, wide containment gains ≥ 14 days mean TTP with 95% lower bound > 0.

Both must hold for the claim "robust within the tested simulator". If H1 holds and H2 fails, report "non-inferior but no shown gain under sampled uncertainty". If H1 fails, name the failing suites.

## Materials

- Python simulator `src/trackc_sim.py` (RK4, dt 0.25 day, exact on/off switch, bisection event times; validated vs SciPy DOP853 within 0.0006 days in script 159).
- Known-parameter oracle `src/trackc_oracle.py` (illustration only).
- No patient data. No downloads. A laptop (8 logical CPUs). Script 168 needed 6,648 s for 17 suites × 5 arms × n 150 (`sealed_results.json` field `runtime_s`).

## Procedures

1. **Freeze arms (no retuning).** Wide containment (p* 0.10, horizon 2 cycles, m prior [3, 20]) and narrow containment (p* 0.05, h 2) from `output/trackc_sealed_containment/frozen_protocol.json`. Comparators: raw single threshold 0.6 (pre-set primary comparator), raw on/off 0.6/0.5, front-loaded periodic, Zhang-style rule (1.0/0.5). Ablation: filtered-burden threshold without the risk rule, same wide prior, tuned on the script-167 tuning draw only, before any new scoring.
2. **New misspecification axes (not used in 150–168).**
   - Acquired resistance: sensitive cells convert to resistant at rate μ ∈ [1e-4, 1e-3] per day under drug.
   - Drug-strength drift: m falls by 2–5% per treated course.
   - Missed readings: each reading skipped with probability 0.2 / 0.4 (controller carries forward).
   - Delayed readings: reading reflects burden 7 days before the decision.
   - Different endpoint: progression at 1.25·B0 instead of 1.2·B0 (controller unchanged).
   - Out-of-range start burden: B0 ∈ [0.5, 0.7].
   Each at two pre-set severities → 12 suites, plus a nominal reference.
3. **Sampled-uncertainty suite (H2).** For each virtual patient, draw independently: noise multiplier ∈ {1, 2, 3}, bias ∈ {0, +10%, +20%}, drug strength m from [3, 20] log-uniform, growth r from [0.001, 0.10] log-uniform, competition (1,1) or (2, 0.5), and one new axis from step 2 with probability 0.5. Weights fixed now; not tuned.
4. **Size.** n 600 scenarios per suite (4× script 168) and 1,500 for the sampled suite, to narrow intervals in noisy suites. Two noise replicates per scenario. Rough runtime estimate from 168: about 4–6 hours on this laptop; to be measured, not assumed.
5. **Seeds.** New scenario seeds 20261300 + k and noise seeds 1_300_000 + 10k + {1, 2}; never used before. Filter Monte Carlo: three particle-filter seeds per scenario, averaged, to reduce resampling sensitivity.
6. **Seal.** Write the protocol in code with a SHA-256 hash, commit and push it before the first scoring run, and refuse to overwrite the result file. A smoke run at reduced n writes elsewhere and is not interpreted.
7. **Score once.** Report every suite, every arm, with intervals, whether or not it helps the claim.

## Risk Assessment

Computer simulation only. No human subjects, animals, patient data, or wet lab. Risk is misinterpretation: the output must not be read as a treatment recommendation, a drug result, or evidence of clinical benefit.

## Data Analysis

- Primary statistic: mean capped TTP difference vs raw single 0.6, per suite. Two-stage paired bootstrap (resample scenarios, then replicates), 4,000 draws, seed fixed in protocol.
- H1: all 13 suites must have lower 95% bound > −7 days. No multiplicity adjustment is needed for an all-suites rule (it is already conjunctive and conservative).
- H2: sampled suite lower 95% bound > 0 and point gain ≥ 14 days.
- Secondary (descriptive, no claims): CVaR_0.10 by both estimands (`cvar_scenario_mean`, `cvar_joint`), courses used, capped fraction, comparison vs front-loaded and vs the filtered-threshold ablation, narrow vs wide trade-off.
- The oracle may be shown as an upper-bound illustration only.

## Limits that will be stated with any result

1. One assumed ODE family; all parameter ranges are assumptions, not measurements.
2. The misspecification axes are still chosen by the author; the sampled suite reduces but does not remove this.
3. The rule's hypothesis came from post hoc exploration (scripts 150–158); 167–169 are sealed re-tests, not an independent replication by another group.
4. The comparator set is small; a better non-filter rule may exist.
5. Front-loaded dosing has better lower-tail outcomes in several settings.
6. Filter latent-resistance estimates are overconfident (script 166).
7. Simulated TTP is not RANO progression or survival.
8. CFB-GBM or other imaging cohorts could test only a separate retrospective forecasting question; they cannot validate this policy or any counterfactual treatment benefit.

## Bibliography

West, J., Adler, F., Gallaher, J., Strobl, M., Brady-Nicholls, R., Brown, J., Roberson-Tessi, M., Kim, E., Noble, R., Viossat, Y., Basanta, D., & Anderson, A. R. A. (2023). A survey of open questions in adaptive therapy: Bridging mathematics and clinical translation. *eLife, 12*, e84263. https://doi.org/10.7554/eLife.84263

Liu, Z., Zhang, J., Hong, L., Nie, Q., & Sun, X. (2025). Multiscale mathematical model-informed reinforcement learning optimizes combination treatment scheduling in glioblastoma evolution. *Science Advances, 11*(32), eadv3316. https://doi.org/10.1126/sciadv.adv3316

Zhang, J., Cunningham, J. J., Brown, J. S., & Gatenby, R. A. (2017). Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer. *Nature Communications, 8*, 1816. https://doi.org/10.1038/s41467-017-01968-5

Liu–West particle filter and prior-art references for belief-based dosing: to be added after the prior-art search.
