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
