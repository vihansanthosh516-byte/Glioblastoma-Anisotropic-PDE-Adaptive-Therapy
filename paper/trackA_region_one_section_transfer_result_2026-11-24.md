# Track A region: one labelled section to other sites — RESULT (2026-11-24)

One locked run (after Deviation 1, a memory-only fix; see the locked plan). Code SHA-256 1034dbe1...6ef3, matching commit f5cc479. Output `output/trackA_region_transfer/result.json` (44 source-target pairs; target sections 1,082 to ~4,000 scored spots, median 2,710).

## Frozen rule: FAIL

Two reasons:
1. Only **3** donors were eligible (AT10, AT3, AT4). AT13 did not have >= 20 core and >= 20 edge spots at two sites. The rule needed at least 4.
2. In AT3 the candidate lost to its ablation (0.9635 vs 0.9703).

| Donor (pairs) | Candidate (rank-normalized) | Ablation (no rank) | kNN transfer | Saadawy markers | Bastola markers | Best niche | Cell-state model |
|---|---:|---:|---:|---:|---:|---:|---:|
| AT10 (4) | **0.970** | 0.944 | 0.910 | 0.408 | 0.718 | 0.712 | 0.857 |
| AT3 (4) | 0.964 | **0.970** | 0.950 | 0.267 | 0.854 | 0.868 | 0.936 |
| AT4 (36) | **0.897** | 0.882 | 0.837 | 0.443 | 0.839 | 0.714 | 0.671 |

Donor-mean AUROC, core vs edge on sections from other sites.

## What the numbers show (descriptive; the rule failed)

- **One labelled section transfers well to other sites of the same tumor.** Both expression classifiers reached AUROC 0.88-0.97 on distant sections in all three donors. This is the main positive observation, but it was not the pre-registered claim (which required the rank-normalized candidate to win everywhere).
- The new step (rank normalization) helped in 2 of 3 donors by 0.015-0.026 and lost in 1 by 0.007. Its added value is not shown.
- Trained per-patient models beat untrained published marker scores by a wide margin. The Saadawy leading-edge vs core signature scored below 0.5 in every donor (0.27-0.44), i.e., it ranked in the opposite direction on these labels.
- **Secondary: confidence did not transfer.** Conformal sets built for 90% coverage reached only 59% coverage on other sites when calibrated on the source site's replicate section (44% with in-sample calibration). Calibration on one site does not hold at another site of the same tumor.

## Weakest points

- 3 donors; AT4 contributes 36 of 44 pairs. Within-patient only.
- PASTE2 alignment comparator not run.
- High AUROCs may partly reflect that core and edge differ strongly in cell composition; the 50% overlap thresholds drop mixed spots, which makes the task easier than real margins.
- Labels and data are now spent for this claim. No retuning.

## Possible next questions (new pre-registration needed, new data needed)

- Why do calibrated confidence sets break across sites (site shift in label balance or composition)? A method that fixes site-level coverage would be a distinct, testable contribution, but it needs fresh labelled sections (e.g., AT15 with a site parser, or a new cohort).
- The opposite-direction result for a published edge signature on pathologist labels is worth a careful check before reporting.
