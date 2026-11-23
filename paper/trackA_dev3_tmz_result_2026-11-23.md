# Track A TMZ development falsification 3: result

2026-11-23. Generated data only. Plan and code frozen in commit ad34b03 before the run ([plan](trackA_dev3_tmz_plan_2026-11-23.md)). One run. Output `output/trackA_dev3_tmz/result.json`, SHA-256 `70c2fead6e76a487e91d6fed437746675f2a7b387ee8bf0165502fe53e37756a`. 20 replicates; replicate-level 95% intervals.

| Scenario | Base ridge | Candidate | Unchecked direction | Shuffled links | Outcome-blind | One-slope contrast | Admitted | Base − cand. CI | Shuffled − cand. CI |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| Confounded transfer | 0.4153 | 0.4077 | 0.3995 | 0.4169 | 0.4137 | **0.4018** | 0.50 | 0.0027 to 0.0126 | 0.0018 to 0.0167 |
| Confounded both | 0.6402 | 0.6326 | 0.6293 | **0.6264** | 0.6397 | 0.6400 | 0.55 | 0.0008 to 0.0146 | −0.0100 to −0.0024 |
| Strong confound | 0.6936 | 0.6923 | 0.6853 | 0.7558 | 0.6940 | **0.6645** | 0.45 | −0.0035 to 0.0062 | 0.0467 to 0.0804 |
| Plain transfer | 0.3334 | 0.3194 | 0.3115 | 0.3115 | 0.3315 | 0.3279 | 0.60 | 0.0082 to 0.0198 | −0.0130 to −0.0027 |
| Null | 0.3117 | 0.3117 | 0.3118 | 0.3118 | 0.3121 | 0.3113 | 0.10 | 0 to 0 | 0 to 0.0001 |

Origin MAE, lower is better. Bold = best arm (unchecked direction ties or is close in some rows).

Rules: confounded transfer beats base — pass; beats shuffled — pass. Strong confound beats base — **fail**; beats shuffled — pass. Confounded both safe — pass. Null safe — pass. Overall: **fail** under the frozen rule.

## What this shows

1. **True sister links now carry information that shuffled links do not.** Under training-only confounding, the shuffled-link direction is worse than the true-link candidate (strong confound: 0.0467 to 0.0804 MAE). This is the first fixture where link destruction hurts. It supports the mechanism in principle, in generated data only.
2. **The candidate's replication gate is too strict.** It admitted the direction in only 45–60% of draws with real signal. The unchecked direction beat the candidate in every signal scenario here. In run 2 the same gate was what protected the reversal case. So the gate trades gain for safety, and this tradeoff is not yet tuned or justified.
3. **The old one-slope contrast is the best arm under confounding** (0.4018 and 0.6645). It also uses true links. So the specific "replicated direction" form is not the best way to use sister links in this fixture.
4. When the confound also exists in new origins, the confound becomes real signal, and shuffled links win (0.6264). Expected.

Weakest point: one hand-built confound structure. Real origin-level confounding in TMZ assays is not measured. No biology or power claim follows.

## Next development step (not run)

Freeze one combined candidate: one-slope contrast weight chosen by a replication check across disjoint origin halves, falling back to origin-balanced ridge when the check fails. Test on run-2 reversal plus run-3 confounding fixtures together, with fresh seeds. Any real gain claim still needs origins with at least two independently derived, assay-matched cultures (Park / Moniz Garcia counts pending).
