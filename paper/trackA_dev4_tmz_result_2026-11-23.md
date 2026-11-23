# Track A TMZ development falsification 4: result

2026-11-23. Generated data only. Plan and code frozen in commit cb38cdf before the run ([plan](trackA_dev4_tmz_plan_2026-11-23.md)). Fresh seed 202611234. One run. Output `output/trackA_dev4_tmz/result.json`, SHA-256 `df837f976517e4ef1d2f182e0e162f646a67ef3bae9e88c5ce39ee7d644dae3f`. 20 replicates; replicate-level 95% intervals.

| Scenario | Base ridge | **Gated contrast** | Contrast unchecked | Shuffled-link contrast | Replicated direction | Admitted | Base − cand. CI | Shuffled − cand. CI |
|---|---:|---:|---:|---:|---:|---:|---|---|
| Reversal | 0.9165 | 0.9198 | 0.9240 | 0.9130 | 0.9165 | 0.30 | −0.0066 to −0.0001 | −0.0209 to 0.0074 |
| Confounded transfer | 0.4564 | 0.4465 | 0.4402 | 0.4912 | 0.4479 | 0.65 | 0.0062 to 0.0135 | 0.0328 to 0.0565 |
| Strong confound | 0.6262 | 0.6065 | 0.5980 | 0.7229 | 0.6189 | 0.70 | 0.0127 to 0.0267 | 0.0931 to 0.1397 |
| Confounded both | 0.6376 | 0.6370 | 0.6358 | 0.6356 | 0.6296 | 0.60 | −0.0006 to 0.0018 | −0.0106 to 0.0078 |
| Plain transfer | 0.3525 | 0.3444 | 0.3444 | 0.3487 | 0.3410 | 1.00 | 0.0061 to 0.0102 | −0.0019 to 0.0105 |
| Unequal sisters + confound | 0.4187 | 0.4150 | 0.4109 | 0.4726 | 0.4141 | 0.45 | 0.0016 to 0.0058 | 0.0457 to 0.0694 |
| Null | 0.3009 | 0.3008 | 0.3007 | 0.3062 | 0.3009 | 0.75 | −0.0004 to 0.0006 | 0.0020 to 0.0088 |

Origin MAE, lower is better.

**Frozen rule: PASS, 9 of 9** (beats base and shuffled links in all three confounded-transfer scenarios; within −0.01 of base in reversal, confounded-both and null).

## What this means

- In generated data with origin-level confounding that does not carry to new origins, true sister links give a real gain. Destroying the links costs 0.03 to 0.14 MAE. This is the first development result where the sister-link mechanism is identifiable against both comparators.
- This is a **development pass only**. It does not touch the Positive, Unique or Substantial boxes, which need real eligible data.

## Limits (state these with the result)

1. **Post hoc gate.** The "sum of half gains" rule was chosen after run 3 showed the stricter rule admitted too rarely. Run 4 used fresh seeds, but the rule itself is not pre-specified from first principles.
2. **The gate costs gain.** The unchecked contrast beat the gated one in every confounded scenario (e.g. 0.5980 vs 0.6065). The gate's benefit is limited: in reversal it cut the loss versus base from 0.0075 to 0.0033, but the candidate is still slightly worse than base there (interval −0.0066 to −0.0001).
3. **Hand-built confound.** Whether real TMZ assays have origin-level confounding that does not transfer is unknown. If it is absent (plain transfer), the gain over base is small and shuffled links tie.
4. **Close prior art.** The one-slope contrast reduces to origin-balanced ridge at one contrast weight (run 1). The gain here comes from a different weight plus true links. Within/between-group regression (Mundlak-type) is established. The novel part to test is the replication-gated use of independently derived sister cultures for drug-response prediction on new origins. A full prior-art search for this exact use is not yet done.
5. One fixture family, 40 training origins. No biology, power or clinical claim.

## Next dependency

Real origins with at least two independently derived, assay-matched cultures with untreated RNA and measured TMZ viability (Park / Moniz Garcia counts pending). Then: full prior-art search, power estimate for the actual rule, and the signed external lock.
