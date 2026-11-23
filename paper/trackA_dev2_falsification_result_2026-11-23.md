# Track A development falsification 2: result

2026-11-23. Generated data only. Plan frozen first ([plan](trackA_dev2_falsification_plan_2026-11-23.md), commit 9a77af6). Code commit 37c565f, `src/trackA_dev2.py`. One run. Output `output/trackA_dev2/result.json`, SHA-256 `85ea206bb4745c686c858e53acbef3d1d17e384fc9fdc897cf163c1dc834fdad`. 20 replicates per scenario, 160 evaluation pairs or 100 evaluation origins each. Intervals are replicate-level normal intervals. Earlier run outputs are kept unchanged.

## Starting-point check

- The implementation files from run 1 were untracked in Git. They are now committed (d4adbda) without changes.
- Scratch copy: 13 of 13 development tests passed. The benchmark output reproduced SHA-256 `DDE30990...442C1` byte for byte.
- Gmail of the author shows no replies from GBM-Space, Petti, Park or Moniz Garcia teams. The sent requests are not visible in this mailbox, so they may have gone from another account. Status of the four requests: no reply on record.

## Regional: harm-calibrated override — FAILED

| Scenario | Candidate | Unchecked override | Nonspatial | Strongest other (arm) | Candidate − strongest, 95% CI | Overrides / harmed / improved |
|---|---:|---:|---:|---|---|---|
| Clean | 0.924 | 0.924 | 0.727 | 0.925 (shared nominal) | −0.003 to 0.003 | 777 / 73 / 704 |
| Shared error | 0.709 | 0.234 | 0.709 | 0.796 (molecule) | −0.101 to −0.074 | 0 / 0 / 0 |
| Divergent error | 0.718 | 0.733 | 0.718 | 0.999 (alternative) | −0.295 to −0.267 | 0 / 0 / 0 |
| Calibration mismatch | 0.238 | 0.231 | 0.705 | 0.797 (molecule) | −0.580 to −0.537 | 1588 / 1541 / 47 |
| Density only | 0.713 | 0.713 | 0.713 | 0.715 (molecule) | −0.007 to 0.003 | 0 / 0 / 0 |
| Missed contamination | 0.702 | 0.708 | 0.713 | 0.801 (molecule) | −0.114 to −0.084 | 938 / 487 / 451 |

Rules: clean gain over strongest — fail (tie). Harm ≤ 0.02 per patient in every adverse case — fail (0.48 in calibration mismatch, 0.15 in missed contamination). Safeguard identifiable — pass (shared error: candidate 0.709 vs unchecked 0.234; the check rejected in 20 of 20 draws).

Diagnosis. The harm check works when the calibration data contain the same error as the test. It cannot see an error that is absent from calibration. That is the same failure as run 1. A threshold rule does not fix it. The limit is measurement: the reference must observe the shared error in the evaluation condition. When the check passes in clean data, the candidate only equals the plain shared spatial model. No spatial contribution claim is earned.

## TMZ: replicated within-origin direction — MECHANISM UNPROVEN

| Scenario | Base ridge | Candidate | Unchecked direction | Shuffled-link direction | Outcome-blind | One-slope contrast | Admitted | Base − cand. CI | Shuffled − cand. CI |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| Transfer | 0.3428 | 0.3224 | 0.3160 | 0.3170 | 0.3419 | 0.3366 | 0.75 | 0.0129 to 0.0279 | −0.0111 to 0.0003 |
| Reversal | 0.8756 | 0.8756 | 0.8877 | 0.8874 | 0.8771 | 0.8836 | 0.00 | 0 to 0 | 0.0067 to 0.0170 |
| Unequal sisters | 0.3109 | 0.2945 | 0.2920 | 0.2927 | 0.3092 | 0.3075 | 0.85 | 0.0112 to 0.0215 | −0.0049 to 0.0013 |
| Assay shift | 0.3313 | 0.3110 | 0.3088 | 0.3091 | 0.3294 | 0.3266 | 0.85 | 0.0122 to 0.0285 | −0.0045 to 0.0008 |
| Null | 0.3148 | 0.3148 | 0.3148 | 0.3148 | 0.3150 | 0.3136 | 0.00 | 0 to 0 | 0 to 0.0001 |

Origin MAE, lower is better. Rules: beats base in transfer — pass. Beats shuffled links in transfer — fail. Safe in reversal and null — pass.

Diagnosis. The candidate beats origin-balanced ridge and the old one-slope contrast. But the shuffled-link direction does as well when between- and within-origin effects agree. In that case any direction estimate finds the same genes, so true sister links add nothing to accuracy. True links matter only in the reversal case. There the replication gate rejected the direction in 20 of 20 draws and kept base accuracy, while the shuffled-link direction lost 0.0067 to 0.0170 MAE. So the measured value of true links in this fixture is safety under reversal, not gain. A gain claim needs a design where the within-origin direction transfers but is hard to find without true links. That is a fixture design question, not yet a result. Weakest point: one generated fixture family; no biological effect or power can come from it.

## Boxes (both claims)

| Claim | Positive | Unique | Substantial |
|---|---|---|---|
| Regional | Blocked: no eligible data. Development candidate failed. | Fail in development: no gain over shared spatial model. | Blocked: no reference that sees shared error in the test condition. |
| TMZ | Blocked: no eligible paired origins. | Not established: ties shuffled links on gain; gate shows safety value only. | Blocked: no independent assay or lab. |
