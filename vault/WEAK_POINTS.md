# Weak-point register (fix plan)

Created 2026-10-05. Rule: every weak point gets a fix, a phase, and a status. Do not delete rows. Status: Open, In progress, Fixed (rerun shows it), Accepted (cannot fix; stated in paper).
Numbers: [[PAPER_FINDINGS_LEDGER]]. Limits: [[LIMITS_CURRENT]].

| ID | Weak point | Evidence | Fix | Phase | Status |
|---|---|---|---|---|---|
| W1 | Primary population mixes WHO grades (30 of 134 not GBM) | consort_mu.json | GBM-only sensitivity for every endpoint (Amendment 2) | 2-3 | Rule fixed; run pending |
| W2 | Informative follow-up: eligible patients progressed 88.8% vs 47.8% | cohort_tables.json | Stratify by progression; inverse-eligibility weighting; wording "among patients with usable pairs" | 3 | Rule fixed; run pending |
| W3 | Parameter fit sits on grid edge (4/5 folds rho 0.1, d 0.01) | forecast_labels_core/results.json | Extended grid on development data, before freeze | 2 | Done 2026-10-07: PDE re-run finished; result negative (ledger 2w) |
| W4 | Median Dice gain only +0.005 | results.json | Patient-level distribution, catastrophic-error rate, subgroup gain, selector | 3-4 | Open |
| W5 | Gain comes from growth, not direction (iso_homog +0.028 vs aniso +0.026) | results.json | Ablation ladder, volume-matched baseline, geometric baseline | 2 | Open |
| W6 | Script 81 folds differ from manifest folds | script 81 vs 94 | Rerun all forecasts on `data/manifests/split_mu.csv` | 2-3 | Fixed (script 100 used the manifest split) |
| W7 | Oracle baseline (`uniform_dilation` uses true scan-2 volume) | run_improved_aniso.py | Replace with causal geometric baseline; keep oracle labelled as upper bound | 2 | Fixed (causal geometric baseline in script 98; oracle labelled upper bound) |
| W8 | Core size checked at native resolution, plan says 2 mm | script 94 | Recompute at 2 mm in the forecast loader; report patients whose status changes | 2 | Fixed (6 of 334 pairs have an empty core at 2 mm: excluded from ladder and PDE, listed) |
| W9 | Track A has 21 patients; folds hold 3-5; SD up to 13 points | cgat_leak_free.json | Leave-one-patient-out with patient-bootstrap CI; pseudo-bulk analysis; look for a second public GBM single-cell cohort with region labels | 6 | Open |
| W10 | "Healthy" = 3 non-GBM donors, all Fresh | trackA_audit.json | Rename; drop Healthy from claims; test tumour-region contrast only within GBM patients (Core vs Periphery, 8 patients with both) | 6 | Open |
| W11 | cells-per-patient alone scores 64.6% (shortcut) | trackA_audit.json | Rerun baselines on a cell-count-balanced subsample (equal cells per patient) | 6 | Open |
| W12 | Graph and cVAE built using all cells and labels | trackA_audit.json | Use only gat_pca_clean-style pipelines; rebuild cVAE on training patients only if kept | 6 | Open |
| W13 | Kill rates assumed, same for all 21 real-growth patients | ledger C1 | Latent per-patient k with prior; posterior update; 2x2 heterogeneity design | 5 | Open |
| W14 | D not identifiable (124/154 at lower bound) | LIMITS L5 | Synthetic recovery, profile likelihood, report identifiable combination | 4 | Open |
| W15 | LUMIERE masks not downloaded (30 GB); external test cannot run | LIMITS L1 | Download masks before the freeze (needs approval and disk space); until then external claims are volume-only | 3 | Open |
| W16 | LUMIERE volumes noisy (median no-change error 0.8-1.8 log units) | ledger 2l | Segmentation-perturbation noise floor; abstain when SNR low | 3-4 | Open |
| W17 | rho > 0.02 cutoff chosen after data, n=10 | LIMITS L16 | Continuous growth-rate analysis; train-locked threshold with bootstrap CI | 5 | Open |
| W18 | MGMT null is underpowered (needs ~2,807 events) | limit_checks.json | State minimum detectable HR; use MGMT as context only | 6 | Accepted |
| W19 | Zone-expression files for scripts 43-44 likely synthetic | LIMITS L32 | Do not use; if Track A to B link is wanted, rebuild from real data with patient IDs that match | 6 | Open |
| W20 | No prospective validation | LIMITS L22 | Not possible. Wording: retrospective, in silico | all | Accepted |
| W21 | Clamp u >= 0 adds mass under strong off-diagonal anisotropy (2.8% in a stress case) and is not logged | solver_verification.json V5 | Log clamp mass per forecast; report per-patient clamp gain; rerun aniso arm at h = 1 mm for a subset | 2-3 | In progress (clamp mass logged per run in script 100) |
| W22 | At h = 2 mm the asymptotic front speed is off by up to +50% in the cell most folds chose (D 0.01, rho 0.1); radius error 0.64 mm at 70 d | solver_verification.json V7b, V9 | Extended-grid check (W3) run at h = 2 mm and h = 1 mm on a 30-patient development subset; report the shift | 2 | Done 2026-10-08: pilot shift +0.0047 at reference cell, best cell unchanged; A6.3 rule not triggered, 2 mm stays primary (ledger 2z) |
| W23 | The PDE gain (+0.026) may be partly matched by a no-PDE population-rate growth baseline (+0.0136; first pairs +0.0159; corrected from +0.031 which was wrong) | baseline_ladder.json | Rerun PDE arms on manifest folds, same pairs; paired comparison PDE vs geometric_train_rate is the real test of the PDE | 3 | Answered 2026-10-07: PDE vs geometric -0.0019 (GBM-only), not better (ledger 2w) |
| W24 | Noise floor in the SNR analysis depends on voxel size (2 mm) | forecastability.json | Recompute at 1 mm; add a measured floor from segmenter disagreement (LUMIERE has two segmenters) | 3 | Open |
| W25 | Baseline bootstrap used 2,000 draws; plan says 10,000 | baseline_ladder.json | Final primary run uses 10,000 | 3 | Fixed for scripts 98 and 100 (ledger 2y); script 102 rerun running |
| W26 | 5 older tests fail (test_hybrid_controller x?, test_inverse_estimation): module `vcs` has no `GbmTherapyEnv`; inverse-estimation KeyError. Not touched by Phase 1-2 changes (new files only) | pytest tests | Triage each: fix or retire; CLAUDE.md claims 56 passing | 3 | Fixed 2026-10-05 (93 tests pass; see ledger 2u) |
| W27 | Manifest flags PatientID_0007 ineligible (one scan has no day) but its 2 valid pairs are analysed; 2 eligible patients with an empty-core pair are dropped; analysed n = 133 (105 GBM), not 134 (104) | ledger 2x; oof_pairs.csv | Decide in the open: keep the pair-table rule (documented). Do not edit the outputs | 3 | Rule fixed in Amendment 6 A6.1; CONSORT figure to show both counts |
