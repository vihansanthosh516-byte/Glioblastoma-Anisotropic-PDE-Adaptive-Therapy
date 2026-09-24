# Project State — Last updated 2026-09-23

## What this project is
Computational glioblastoma modeling with three tracks. Focus is currently 
on Track B (adaptive therapy).

## Verified results (committed)
- **Script 27** (a2802ac): ABA lattice, front velocity 0.118-0.120, sustained
- **Script 28** (7bd7def): FK-PDE, numerical 3.62 vs analytical 4.00, 9.4% error
- **Script 29**: integrated invasion, free-propagation 26.6 um/hr vs analytical 28.3
- **Script 43** (e76d26a): stromal feedback, per-patient rho, tumor mass 14x variation
- **Script 44** (17bfade): adaptive therapy, 61 real patients. 
  Resistance: median 99.2% MTD vs 30.6% adaptive. 
  54/61 select more resistance under MTD, 0/61 under adaptive.

## Currently running
- Threshold test: rerun script 44 with THRESHOLD_OFF=0.50 to see if high-rho 
  failure is a setpoint artifact or real biology

## Open questions
- Is high-rho failure (10/61 progressed earlier under adaptive) real biology?
- Where do the 42 negative-rho responders go in the paper?
- Scripts 46, 47, 48 still need per-patient param patches

## Constraints
- E_MAX_RATIO=1000 (empirical sweep: 500x undershoots, 2000x saturates)
- THRESHOLD_OFF=0.80 (under test)
- Negative-rho patients excluded from adaptive therapy sim
- All scripts use PROJECT_ROOT/OUTPUT_DIR pattern

## Where things live
- src/44_adaptive_therapy.py
- output/adaptive_cohort_summary.json (headline)
- output/adaptive_geometry_metrics.json (61 per-patient records)
- TRACK_A_RESULTS.md (scratchpad)