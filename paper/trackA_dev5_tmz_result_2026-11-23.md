# Track A TMZ development falsification 5: result

2026-11-23. Generated data only. Plan and code frozen in commit a528f5f before the run ([plan](trackA_dev5_tmz_plan_2026-11-23.md)). Fresh seed 202611235. One run. Output `output/trackA_dev5_tmz/result.json`, SHA-256 `8753a96a6ebff2b6696c15eb4ae32826cb94d0e20b849f578cbcb5c443d95b33`. 20 replicates; replicate-level 95% intervals.

| Scenario | Base ridge | Gated contrast (cand.) | Contrast unchecked | Shuffled-link contrast | Mundlak ridge | **Within-only** | Strongest − cand. CI |
|---|---:|---:|---:|---:|---:|---:|---|
| Reversal | 0.9045 | 0.9077 | 0.9126 | 0.9067 | 0.9160 | **0.8697** | −0.0573 to −0.0187 |
| Confounded transfer | 0.4224 | 0.4112 | 0.4065 | 0.4689 | 0.4841 | **0.3196** | −0.1157 to −0.0673 |
| Strong confound | 0.6213 | 0.6062 | 0.5943 | 0.7192 | 0.7627 | **0.4271** | −0.2359 to −0.1223 |
| Confounded both | 0.6436 | 0.6411 | 0.6406 | 0.6342 | **0.6332** | 0.7419 | −0.0170 to 0.0013 |
| Plain transfer | 0.3521 | 0.3457 | 0.3457 | 0.3463 | 0.3319 | **0.3269** | −0.0310 to −0.0068 |
| Unequal sisters + confound | 0.3956 | 0.3912 | 0.3872 | 0.4382 | 0.4698 | **0.3173** | −0.0927 to −0.0551 |
| Null | 0.3120 | 0.3114 | 0.3115 | 0.3171 | 0.3221 | **0.3006** | −0.0148 to −0.0068 |

Origin MAE, lower is better. Strongest comparator = best of base ridge, Mundlak ridge, within-only, shuffled-link contrast.

**Frozen rule: FAIL, 1 of 6.** The candidate lost to the strongest established comparator in every confounded-transfer scenario, in reversal and in null. It passed only in confounded-both.

## What this means

- The gated contrast's dev4 pass was a pass against weak comparators. An established within-family-style predictor (within-origin coefficients, then a one-number recalibration on origin means) beats it by 0.07 to 0.24 MAE wherever origin-level confounding does not transfer. It also wins in reversal, because it never uses the between-origin direction.
- Mundlak ridge did poorly when the confound does not transfer: the origin-mean features carry the confound to new origins.
- Within-only loses when the confound is real signal in new origins too (confounded both: 0.7419 vs 0.6332). That is its known weakness.
- **The gated contrast is retired as a candidate.** Results are kept.

## Where this leaves the TMZ method

- The best arm is an existing method (family-GWAS direct-effect logic, applied to sister cultures). Applying it to independently derived tumour cultures for drug response may be new as an application (unverified; see [prior art](trackA_tmz_gated_contrast_prior_art_2026-11-23.md)). It is not a new method.
- A genuinely new method would need to beat within-only where it is strong AND Mundlak / base ridge where the confound transfers (confounded both). One option to test next: choose between within-only and pooled ridge per training set by the same disjoint-halves replication check. That is a model-selection rule on established estimators; its novelty would also need checking.
- Weakest point of all TMZ development results: hand-built confound; real TMZ data have not shown whether origin-level confounding exists or transfers.
