# Track C final status report and plan amendment (2026-10-31)

SIMULATION ONLY. All parameters, including switching (mu) and reversion (nu) rates, are assumed. Time to
progression (TTP) here is a capped model quantity, not patient progression or survival. No clinical benefit,
drug efficacy or novelty is claimed. Every number below was read back from the cited JSON on 2026-10-31.

## 1. Verified sealed results (600 scenarios per suite, CVaR10 = mean of the worst 10% of scenario means)

Arms (both scripts): tuned arm `riskf (p*=0.05, horizon 2, front K=4)`; front-loaded `periodic 1`;
comparator `raw 0.6`. Script 173 tuned on development draws and re-selected the 171 arm, so 173 is a
sequential re-test of the same rule on new seeds, not an independent replication of a new rule.

| Script / suite | Arm mean TTP vs front-loaded (d, 95% CI) | H1: CVaR10 arm - front (d, 95% CI; pass if lower > -7) | H2: mean arm - comparator (d, 95% CI) |
|---|---|---|---|
| 171 nominal | +32.5 [28.4, 36.7] | -1.60 [-4.17, -0.50] pass | +22.0 [16.9, 27.4] |
| 171 state_linked_mod | +30.0 [26.0, 34.1] | **-4.18 [-7.97, -1.44] FAIL** | +33.9 [28.1, 40.1] |
| 171 state_linked_sev | +27.8 [23.9, 31.8] | -2.64 [-5.54, -0.68] pass | +54.2 [46.4, 62.5] |
| 171 mixture | +6.3 [1.7, 10.4] | 0.00 [0.00, 0.00] pass (identical tails) | +17.8 [10.2, 25.2] |
| 173 nominal | +30.6 [26.5, 34.6] | **-3.78 [-8.04, -1.33] FAIL** | +17.1 [13.1, 21.3] |
| 173 state_linked_mod | +30.3 [26.5, 34.2] | -1.81 [-3.51, -0.72] pass | +35.9 [30.2, 42.0] |
| 173 state_linked_sev | +29.9 [25.8, 33.9] | -2.23 [-5.98, -0.52] pass | +59.7 [51.8, 68.5] |
| 173 mixture | +8.7 [4.3, 13.0] | 0.00 [0.00, 0.00] pass (identical tails) | +16.0 [9.7, 22.4] |

Sources: `output/trackc_sealed_tail_hybrid/sealed_results.json` (SHA-256 e0e2ecad...), verdict H1 fail
[state_linked_mod], H2 pass, H3 pass; `output/trackc_sealed_tail_trigger/sealed_results.json` (SHA-256
b004ca91...), verdict H1 fail [nominal], H2 pass.

Reading:
- Mean TTP and worst-decile performance point in opposite directions. The arm gained mean days in every suite,
  but its CVaR10 point estimate was at or below front-loaded in every suite, and every non-mixture CI upper
  bound was below zero.
- H1 (lower-tail non-inferiority, -7 d margin) failed in both sealed tests, in different suites.
- H2 passed only against the prespecified `raw 0.6` comparator. The "vs front-loaded" mean column is a
  secondary contrast; it does not offset the tail result.
- Mixture suite H1 = 0 [0, 0] because both arms had the same tail values; this is not evidence of protection.
- 171 and 173 are sequential looks at one rule family; their passes are not two independent confirmations.

## 2. Verified development screens (not sealed; no test seeds)

Script 177, `output/trackc_dev_plasticity_sensitivity/dev_sensitivity_results.json` (SHA-256 e0246ecd...):
36 scenarios, 2 noise replicates, 9 assumed (mu, nu) cells in {0, 0.01, 0.04}/day.
- Forced-holiday policy (`riskh`, max run 3): mean above front-loaded in 8 of 9 cells; CVaR10 below
  front-loaded in 7 of 9 cells. Example mu=nu=0.04: mean 743.5 vs 665.7 d; CVaR10 208.2 vs 361.7 d.
- Weak points: 36 scenarios means each CVaR10 rests on about 3.6 scenarios; no CI for the tail contrast;
  rates are assumed; the run used a parallel harness (recorded in `execution_manifest.json`).

Script 178, `output/trackc_policy_choice_upper_bound/development_bound.json`: reuses the 177 outcomes for a
post-hoc, clairvoyant choice between front-loaded and wide containment per scenario. Mean gain over the better
fixed policy 0.62-11.98 d; tail gain 0 in 7 of 9 cells, 3.26 and 16.16 d in two. This is not a new test and
not an implementable policy.

Script 176 (`output/trackc_dev_plasticity/dev_results.json`) earlier gave the same pattern: forced holidays
lowered CVaR10 in every suite (`paper/trackc_mechanism_review.md`).

## 3. Search for independent rate data (2026-10-31, targeted, not exhaustive)

| Source | What it measures | Constrains nu or mu? | Access / terms |
|---|---|---|---|
| Rabé et al. 2020, Cell Death Dis (PMC6944699) | U251 transient TMZ-tolerant state; rechallenge after 12 weeks off drug | No: one cell line, one washout time | Open article, CC BY 4.0 |
| Banelli et al. (PMC4825557) | Resistance partly reversed after >= 30 days washout | No: one washout duration | Article open; no explicit data licence found |
| Kuo et al. (PMC9101568) | Regrowth after rechallenge and drug removal | No: not sensitivity vs washout time | CC BY 4.0 |
| Lasri et al. 2020, R Soc Open Sci | Stochastic MGMT expression, in vitro rates, one line | Partial for mu-like switching in vitro; not human | Open access |
| LUMIERE (Bern BORIS 175955, Sci Data 2022) | 91 GBM patients, ~7 MRI time points, RANO, partial MGMT | No: no per-patient dose timing or off-drug interval to identify nu | CC BY |
| MU-Glioma-Post, CFB-GBM (already in project, Track B) | Longitudinal MRI masks | No: standard schedules only; no on/off contrast; CFB 0/107 preflight pass | CC BY 4.0 |

Limited conclusion: this search found no public dataset that measures GBM sensitivity at several off-drug
durations or in patients under varied TMZ schedules. That is what would be needed to estimate nu. This does
not prove no such data exists.

## 4. Decision

No independent biological evidence constrains the switching or reversion rates. Without it, any new simulation
only shows what the assumed rates imply. Results from 176-178 already show that the conclusion flips with the
assumed rates (e.g. holiday tail worse in 7 of 9 cells, better in 2). Picking a rate range that favours a
policy would be tuning on seen outcomes. **No new simulation was run, and none is justified now.**

## 5. Plan amendment (future tense)

1. Track C will be reported as an exploratory simulation study of a mean-vs-worst-decile trade-off inside one
   assumed model family.
2. No further sealed Track C test will be opened. H1 will be reported as failed in 171 (state_linked_mod) and
   173 (nominal). No margin, endpoint, suite or scenario range will be changed.
3. The 600-scenario minimum-CVaR primary outcome in `paper/research_plan_complete_project_appendix.md` §5 was not
   run as written; it will be listed as not executed. `paper/research_plan_trackC_draft.md` stays a draft and is
   not approved by this amendment.
4. A new protocol will be considered only if a source gives a measured off-drug recovery time course (several
   washout durations, ideally patient-derived GBM cells), and an outside reviewer agrees the rates bound the
   simulation before any outcome is seen. It would need fresh seeds, predeclared estimands and a power
   justification for the tail contrast.

## 6. Addendum 2026-11-03: record corrections (no result changed)

Source: read-only audit by 14 Haiku 5.5 experts on 2026-11-03, plus Claude re-checks. `[V]` = re-checked by
Claude against the named file on 2026-11-03 (by `src/trackc_doc_check.py` or by hand). `[A]` = taken from an
expert report or an earlier entry and not re-checked here. Old results, margins and endpoints are unchanged.
Only wording and missing context are added. Restart plan: `paper/trackc_restart_plan_2026-11-03.md`.

1. [A] Sealed H1 failed three times, not two. The 169 H1 (wide containment minus raw 0.6, mean TTP, lower
   bound > -7 d in all 13 suites) failed at missed_sev: -4.7 [-16.6, 7.0]
   (`output/trackc_sealed_new_axes/sealed_results.json`). The 171 and 173 H1 (worst decile vs front-loaded)
   failed as shown in section 1. These are different tests. Script 170 was designed after the 169 failure
   (`START_HERE.md`, script-170 entry), so its pass is not an independent confirmation. Section 5, item 2
   lists only 171 and 173.
2. [V] One tuning draw (seed 20261121, n 300, noise keys 1121001 and 1121002) is recorded in the protocols of
   167, 169, 171 and 173. The tuning values for the 167 wide arm and the 171 arm are identical (700.503...),
   and the 171 and 173 tuning values for the chosen arm are identical (699.417...). So 171 and 173 are
   sequential looks at one arm (4 forced courses, then wide-prior containment, p* 0.05, horizon 2), as
   section 1 says. They are not independent replications.
3. [A] The 167 primary arm is narrow-prior containment (p* 0.05, horizon 2). In script 168 that narrow arm
   lost to the comparator when the drug was weaker than its prior: weak_drug_mod -59.3, weak_drug_sev -81.7
   (-98.1 to -65.5), combo_fast_weak -19.8 (`output/trackc_sealed_wrong_models/sealed_results.json`;
   `START_HERE.md`, script-168 entry). The "no harm" statement for 168 applies to the wide arm only. The
   171 and 173 hybrid uses the wide prior and was not an arm in 168, so its weak-drug behaviour is untested.
4. [A] Holiday direction (section 2): forced holidays gave a higher mean in 8 of 9 cells and a lower (worse)
   CVaR10 in 7 of 9 cells. The phrase "cut the worst-10% CVaR" in `paper/trackc_mechanism_review.md` meant
   "lowered the score" (reworded there). Script 176 gave the same direction in all three suites.
5. [V] CVaR10 of scenario means is biased low for feedback rules at few noise repeats. In the 167
   per-scenario file (400 scenarios, 4 repeats), wide containment gives 171.1, 174.5 and 176.0 at 1, 2 and 4
   repeats; raw single 121.1, 125.5 and 132.6; front-loaded 182.2 at every count, because it ignores
   readings (open loop, identical repeats). Scripts 171 and 173 used 2 repeats (`protocol.n_rep = 2`). The
   tail cost they report (1.6 to 4.2 d) is the same size as this bias, so the H1 contrast is not a clean
   estimate of tail harm. The H1 results stay recorded as FAILED.
   Source: `output/trackc_sealed_containment/sealed_ttp_by_scenario_replicate.csv`.
6. [V] Paired view (same file, wide containment minus front-loaded, scenario means over 4 repeats): lower in
   21.0% of scenarios; mean of the worst 10% of paired differences -42.1 d; mean difference +34.2 d.
7. [V] Share of scenarios with no progression by the 1,456 d horizon (same file): 5.9% narrow containment,
   6.1% wide, 3.4% raw single, 2.2% front-loaded. The mean-TTP contrast is a restricted-mean contrast and
   includes this difference.
8. [V] The comparator `raw 0.6` is the lowest value on its tuning grid (0.6 to 1.18): tuning mean TTP 683.1
   at 0.6 and 682.1 at 0.7 (`output/trackc_sealed_containment/frozen_protocol.json`). The optimum may lie at
   or below 0.6. The curve is flat near the edge, so the shift is probably small. This is untested.
9. [V] The oracle docstring called the oracle an "upper reference" and also said its schedules are
   achievable, a lower bound on the optimum. Corrected on 2026-11-03 in `src/trackc_oracle.py`. The headroom
   numbers (mean 83.5 d; tail 3.4 d, CI 1.26 to 6.93 [A], `output/trackc_headroom/headroom_results.json`) are
   lower bounds on the possible gain, not upper bounds.
10. [V] Per-scenario outcome files exist only for script 167. Scripts 168, 169, 170, 171 and 173 saved summary
    JSON and logs only, so their bootstrap CIs cannot be rebuilt from saved data.
11. [V] The filter rules receive the true per-scenario noise sd (`src/trackc_sim.py` line 14: "sigma known to
    the filter"; `ParticleFilter` is built from `scen["sigma"]` at lines 340 and 375) and, for narrow-prior
    arms, a prior equal to the generator's ranges (line 42). The raw rules receive neither. Nominal sealed
    suites are the best case for the filter rules.
12. [V] Solver mode: the default is `EVENT_MODE = "refine"` (`src/trackc_sim.py` lines 91 to 95) and the 167
    protocol records `"event_mode": "refine"`. Scripts 150 to 158 used "grid" [A] (`START_HERE.md`,
    script-159 entry). The change is not stated in section 1.
13. [V] Commit pins: the commits named in `START_HERE.md` (8762ff0, f2d16bf, 313017f) are freeze and result
    commits. The `git_commit_at_run` field in each JSON is the HEAD at run time (for example d55dcb2,
    60234ab, 2f7c373). They are different commits, not a conflict.
14. [V] Scale: the 167 per-scenario arm means run from 596.0 to 701.9 d. The model uses 28 d scan spacing, 12
    courses, a 1,456 d horizon and progression at 1.2 times the start burden (`src/trackc_sim.py` lines 31
    to 36). Real median progression-free survival is about 7 months (not yet verified; Stage 1 target).
15. [V] Superseded (banner added 2026-11-03): `Claude outputs/trackC_results_report_claude.md` (2026-10-25,
    "first positive result", scripts 150 to 156, no sealed results; the file is untracked) and
    `paper/research_plan_trackC_draft.md` (a draft that reports 169 and 170).
16. [V] Seed registry (`src/trackc_seed_registry.py`, `output/trackc_seed_registry.json`). Script 176 says its
    seeds were "unused before". That holds for its scenario seeds (20262100 and up). Its noise keys
    `1_300_000 + 10k + rep` for k = 0 to 2 equal the keys sealed script 169 used for k = 0 to 2, on
    different scenarios. No sealed outcome depends on this, and no result changed. Seeds 20261800 to
    20262099 (dev) and 20262300 to 20262499 (sealed) were unused by any Track C script on 2026-11-03 and are
    reserved in the registry for the restart plan.
