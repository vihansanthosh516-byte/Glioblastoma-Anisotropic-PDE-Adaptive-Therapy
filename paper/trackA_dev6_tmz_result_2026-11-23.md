# Track A TMZ development falsification 6: result

2026-11-23. Generated data only. Plan and code frozen in commit 5656d4b before the run ([plan](trackA_dev6_tmz_plan_2026-11-23.md)). Fresh seed 202611236. One run. Output `output/trackA_dev6_tmz/result.json`, SHA-256 `88fe4baaec2f5947fbdd6bd5f7e33618a5a126065de2fbdf29ad55df85078097`. 20 replicates.

| Scenario | Within-only | Base ridge | Mundlak ridge | Shuffled-link contrast | Halves selector (cand.) | Picks: within / ridge / Mundlak |
|---|---:|---:|---:|---:|---:|---|
| Reversal | **0.8801** | 0.9217 | 0.9444 | 0.9260 | 0.8970 | 15 / 2 / 3 |
| Confounded transfer | **0.3566** | 0.4604 | 0.5130 | 0.5033 | 0.4492 | 6 / 9 / 5 |
| Strong confound | **0.4586** | 0.6492 | 0.8336 | 0.7487 | 0.7213 | 1 / 11 / 8 |
| Confounded both | 0.7052 | **0.6299** | 0.6302 | 0.6314 | 0.6483 | 2 / 14 / 4 |
| Plain transfer | **0.3173** | 0.3455 | 0.3308 | 0.3447 | 0.3231 | 10 / 6 / 4 |
| Unequal sisters + confound | **0.3384** | 0.4160 | 0.5020 | 0.4658 | 0.4806 | 1 / 5 / 14 |
| Null | **0.3011** | 0.3108 | 0.3234 | 0.3161 | 0.3059 | 14 / 5 / 1 |

Origin MAE, lower is better.

Maximum regret over the 7 scenarios: within-only 0.075, base ridge 0.191, selector 0.263, shuffled-link contrast 0.290, Mundlak 0.375.

**Frozen rule: FAIL.** The selector's worst case (0.263) is worse than always using within-only (0.075).

## Diagnosis

- The selector picks badly exactly where it matters. When the confound is in training only, pooled ridge and Mundlak look better on held-out *training* origins because those origins share the confound. So the selector picks them, and they fail on new origins. In strong confound it picked within-only in 1 of 20 draws.
- This is the identifiability limit written into the plan: training-only data cannot reveal whether origin-level confounding will carry to new origins. No rule that uses training data alone can fix this.
- Always using the established within-only estimator was the safest choice in these fixtures.

## Conclusion for the TMZ method line (development)

Six development runs. The best and safest arm is an established within-family-style estimator applied to sister cultures. The project-built candidates (one-slope contrast, replicated direction, gated contrast, halves selector) all failed against it or against shuffled links. No new TMZ method is supported. The remaining open question is empirical: does real sister-culture TMZ data show origin-level confounding that does not transfer? That needs data from Park, Moniz Garcia, Nam 2023 or Salatino 2026. Choosing between within-only and pooled models needs an external signal about transfer (for example, a second lab or assay), not a training-only rule.
