# Track A region: AT15 fresh test of target-adaptive confidence — RESULT (2026-11-24)

One locked run (plan `trackA_region_at15_prereg_LOCKED_2026-11-24.md`, commit 40a1bb8; code hashes in the result match). Output `output/trackA_region_at15_test/result.json`.

## Frozen rule: FAIL

8 of 9 AT15 sections had labels; 4 had >= 20 core and >= 20 edge spots (A2, B2, C2, G1), giving 12 cross-block pairs.

| Arm | Mean accuracy among confident spots |
|---|---:|
| **Candidate (target mixture + spatial smoothing)** | **0.809** (needed >= 0.90) |
| Matched-coverage source probability | 0.795 |
| Matched-coverage smoothed source probability | 0.800 |
| Ablation: no spatial smoothing | 0.776 |
| Ablation: no mixture (smoothed source probability) | 0.809 |
| Unmatched source probability | 0.786 |

Candidate kept 81% of spots as confident. It beat each matched comparator in 5 of 12 pairs (needed 8).

## What this means

- The development gain (0.903 vs 0.69-0.71 on spent Visium v1 donors) did **not** replicate on a fresh donor and platform. On AT15 the methods are within about 0.03 of each other, and the mixture step added nothing over spatial smoothing alone (0.809 vs 0.809).
- The development result was likely shaped by choosing the design after seeing the spent data, or by differences between Visium v1 and CytAssist v2 sections.

## Weakest points

- One donor, 4 eligible sections, 12 non-independent pairs.
- AT15 is now spent.

## Track A region status

Two locked tests failed (one-section transfer rule; target-adaptive confidence). Robust observation across both: a model trained on one labelled section ranks core vs edge well on other sites of the same tumor. No new method beats simple baselines on fresh data. Further work needs new pathologist-labelled multi-section data.
