# Track A region: fresh test of target-adaptive confidence on AT15 (LOCKED, 2026-11-24)

Locked before any AT15 label value is read (only per-section label presence was checked for four files). Code: `src/trackA_region_at15_test.py`, which imports `trackA_region_transfer.py` and `trackA_region_dev_adaptive.py` (hashes in the result). One run after all nine AT15 files finish downloading: `py -3.14 src/trackA_region_at15_test.py`. Refuses to overwrite `output/trackA_region_at15_test/result.json`.

## Why this test

The locked one-section transfer test (`trackA_region_one_section_transfer_result_2026-11-24.md`) showed that core/edge ranking transfers across sites (AUROC 0.88-0.97) but confidence does not (90% conformal sets covered 59%). Development on the spent donors (`output/trackA_region_dev_adaptive/dev_result.json`, 44 pairs) found that refitting the decision boundary on the target section's own scores, then smoothing over neighbouring spots, raised accuracy among "confident" spots to 0.903 (31/44 pairs >= 0.90), versus 0.690 and 0.713 for source-probability baselines at the same fraction of confident spots. Those numbers were chosen after looking at spent data, so they need a fresh test.

## Data

GBM-Space donor AT15: nine Visium CytAssist v2 sections from tissue blocks A2, B2, C2, C3 (two sections), D3, F1 (two sections), G1. Different platform from the development donors (Visium v1). AT15 was excluded from the earlier tests by the file-name rule and has not been scored. Sections with any label are used; spot labels and the 20-spot class minimum are the same as in the locked test. Site = tissue block code; pairs use a source and a target from different blocks.

## Candidate (frozen)

1. Source model: the locked candidate (within-section rank normalization, PCA 30, balanced logistic regression) trained on one labelled section.
2. Predict every spot of the target section (labels not used).
3. Fit a two-component Gaussian mixture to the target logits; P(core) = posterior of the higher-mean component.
4. Average P(core) over each spot and its 6 nearest spatial neighbours.
5. A spot is "confident" if max(P, 1 - P) >= 0.90; prediction = core if P >= 0.5.

## Comparators

Matched-coverage source probability and matched-coverage spatially smoothed source probability (same fraction of confident spots as the candidate, chosen by distance from 0.5). Also reported: ablation without spatial smoothing, ablation without the mixture, and the unmatched source probability.

## Pass rule

PASS only if all hold: at least 6 cross-block pairs; candidate mean accuracy among confident spots >= 0.90; candidate mean higher than both matched-coverage comparators; candidate better than each matched comparator in at least 2/3 of pairs. If no pair is eligible: INCONCLUSIVE.

## Weak points

- One patient. Pairs share sections and are not independent; the result is a within-patient demonstration on a new platform, not generalization.
- The candidate keeps fewer confident spots (dev: 57% vs 75% for raw source probabilities).
- A two-component mixture assumes both core and edge are present in the target section; a section that is almost all core may break it.
- Generic test-time recalibration and spatial smoothing exist in other fields; novelty for this task was not exhaustively checked.
