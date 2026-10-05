# Script 75 — Probe-then-commit and paced heuristic

**Source:** `output/probe_paced_policies/evaluate.json` (script `src/75_probe_and_paced_policies.py`). Budget 35.0, horizon 365 d, progression factor 1.4. Values read from the JSON on 2026-10-04.

## Paced heuristic (arm `paced_heuristic59`), equal drug budget
| Set | n | Day-90 win | TTP minus Stupp, days (95% CI) | longer / equal / shorter | Wilcoxon p |
|---|---|---|---|---|---|
| cohort64 (synthetic) | 20 | 0% | +48.5 (34.9, 61.7) | 20 / 0 / 0 | 8.8e-05 |
| lhs60 (synthetic) | 30 | 0% | **+52.5** (39.0, 66.1) | 30 / 0 / 0 | 1.6e-06 |
| real_test (real MU-Glioma params) | 21 | 0% | **+2.1** (1.2, 3.3) | 15 / 6 / 0 | 4.2e-04 |
| synth_test (synthetic) | 40 | 15% | +15.7 (7.8, 25.7) | 36 / 4 / 0 | 1.5e-07 |

- Pooled: TTP longer in 101 of 111 patients, shorter in none.
- **Quote the four rows, not "+52.5 d".** The +52.5 d is the lhs60 set only. On real parameters the gain is +2.1 d.
- Paced loses on day-90 volume (0% win on three sets).
- The volume threshold does nothing. The finding is "pacing", not "the heuristic".

## Probe-then-commit (cohort64, L = 1 day probe), day-90 win vs Stupp
| Arm | sigma 0 | sigma 0.10 (3 seeds) |
|---|---|---|
| `probe_rule` | **100%** | 45 / 75 / 65 = **62%** |
| `probe_dagger_cond` | 100% | 40 / 50 / 55 = 48% |
| `probe_dagger_ttp` (TTP-trained) | 70% (TTP +58.9 d) | 60 / 40 / 55 (TTP +46 to +53 d) |
| blind `ppo66` | 25% | 25% |

- Noise-free is 100%. Under realistic noise it falls to about 62%.
- Kill rates are identifiable from daily noise-free volumes, not from clinical volumetry (about 5-10% CV, from `NEGATIVES_REVISITED.md`).
- On lhs60 and real_test the win rate under noise falls to 11-46% (from `NEGATIVES_REVISITED.md`; not re-derived here).
- Frame as an identifiability result, not a deployable policy.

## Label correction (2026-10-04)
`real_test` uses real fitted growth rates but the same assumed kill rates [0, 0.05, 0.08, 0.13] for all 21 patients. Figures now label it "real growth, assumed kill". Median kill/rho is 77.7 there vs 2.7-12.5 in the synthetic sets (script 87), which explains the small +2.1 d gain.
