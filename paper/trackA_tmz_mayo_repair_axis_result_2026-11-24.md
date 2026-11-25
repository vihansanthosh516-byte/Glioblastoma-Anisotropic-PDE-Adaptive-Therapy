# Track A TMZ: repair-axis score on Mayo PDX — RESULT (2026-11-24)

One locked run. Plan: `trackA_tmz_mayo_repair_axis_prereg_LOCKED_2026-11-24.md` (commit 4d0e1ef). Code SHA-256 e3a5300b...320e03 (matches the locked commit). Output `output/trackA_tmz_mayo_transfer/result.json`, SHA-256 74e7e8b2923d611da9f6c63994613ca63a4c72b597405a7d9cf08ee4bf44d570.

## Frozen rule: PASS (narrow)

27 orthotopic PDX (MGMT methylated 11, unmethylated 16). Outcome = log TMZ survival ratio. Partial Spearman adjusting for MGMT status.

| Score | Partial rho |
|---|---:|
| **C2 repair axis (MGMT + MMR expression)** | **0.367** (one-sided permutation p = 0.040) |
| 1. MGMT status | 0.000 (by construction) |
| 2. MGMT expression alone (ablation) | 0.165 |
| 3. Molecular group, best indicator | 0.234 |
| 4. Ntafoulis 102-gene score | 0.112 |
| 5. HGCC transfer ridge, MGMT residual | -0.199 |
| 6. HGCC transfer ridge, raw AUC | -0.234 |
| Random 4-gene swaps (MGMT kept): median / 95th pct | 0.112 / 0.322 |

C2 beat every comparator and 97.5% of random swaps. Raw Spearman (no adjustment) 0.384.

## What this supports

In 27 PDX, a fixed, untrained score combining low MGMT expression with high mismatch-repair gene expression ranked in-vivo TMZ survival benefit beyond MGMT methylation status, better than MGMT expression alone, subtype, a published-signature score and in-vitro-trained transfer models. Adding the MMR term roughly doubled the partial correlation over the ablation (0.165 to 0.367).

## Weakest points (state with the result)

- **Barely significant, small n.** p = 0.040 with 27 models; the 95% random-swap threshold was 0.322, so the margin is small. One or two models could change the verdict. Needs an independent cohort.
- **Not fully outcome-blind.** Median TMZ ratio by MGMT group was seen before the lock (expression was not). MGMT status is adjusted for in the test.
- **Measure changed** from in-vitro viability to in-vivo survival benefit (recorded in the plan, user-approved).
- **Novelty is limited.** MGMT and MMR roles in TMZ response are textbook biology. The result confirms that their expression combines usefully in vivo; it is not a new mechanism. A full prior-art check (e.g., MGMT/MMR expression scores in TMZ-treated cohorts) has not been done, so no "first" claim.
- **In-vitro training failed to transfer.** Both HGCC-trained models scored negative partial rho. That is a negative finding for in-vitro to in-vivo transfer.
- **Zhao 2019 signature not compared** (gene list not obtained).
- Simulation-free, but preclinical (mice). No clinical claim.

## Three boxes

- Positive: narrow pass on one small cohort (not untouched in the strict sense; see above).
- Unique: not established. Ablation shows the MMR term contributes; prior-art check pending.
- Substantial: not tested. Needs a shift test (e.g., the human randomized trial GSE7696: TMZ+RT 52 vs RT 28, not yet used).

Next: prior-art check for MGMT+MMR expression scores; pre-register the same frozen score on GSE7696 (treatment x score interaction on survival) as the independent test. Mayo is now spent for this score.
