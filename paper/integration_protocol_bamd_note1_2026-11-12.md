# BAMD protocol note 1 (2026-11-12, before any survival time is linked to features)

Additions to `paper/integration_protocol_bamd_2026-11-12.md` (SHA-256 0dc8d720...). Decided before Gate 1 results and before any survival analysis.

1. **UPENN survival coding** (from category labels in `clin.csv`; only labels and the 27 censor times were seen, not linked to any feature): event = `Survival_Status == "Deceased"`, time = `Survival_from_surgery_days_UPDATED`. All other statuses ("Alive", "Lost to Follow-up", "Deceased - uncertain date of death") are censored at `Survival_Censor` days. IDH1 "Mutated" excluded; "Wildtype" and "NOS/NEC" kept. KPS is missing for 596 of 671 and MGMT for 348, so neither is in B0 (as the protocol already had); MGMT is used only in a sensitivity model on patients with a known result.
2. **Bootstrap.** The point estimate of each C-index uses 10 repeats of 5-fold CV as written. Each of the 1,000 bootstrap resamples uses one 5-fold split (not 10), because the full design needs about 400,000 model fits per model. CIs may be slightly wider.
3. **Gate 1 permutation null** uses one CV repeat per shuffle while the observed rho averages 10 repeats. This makes the p-value conservative (harder to pass).
4. **Cox fitting**: own Newton-Raphson Breslow Cox in NumPy (tested against lifelines on synthetic data) for speed; ridge 1e-4 for numerical stability in all unpenalized models. LASSO Cox by proximal gradient with the penalty chosen by inner 5-fold CV from 5 values.
