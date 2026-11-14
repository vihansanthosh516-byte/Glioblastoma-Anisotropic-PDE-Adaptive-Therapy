# Integration new method: Evidence-Anchored Reliability Fusion (EARF) — protocol (2026-11-14)

Written before any EARF world has run. No policy outcome has been seen. Simulation only.

## Why
Four real-data MRI->biology tests failed (paper/integration_failure_diagnosis_2026-11-14.md). The Track A molecular->TMZ check also showed no association. The current combined controller (`finite_set_minimax_mpc`) takes the worst case over all 135 models in a fixed synthetic grid, including strong signal loadings (|a|, |b| = 0.5) that the real data do not support. EARF keeps the same controller but lets each signal count only as much as the measured evidence allows.

## Method
1. Measured reliability anchors (fixed now, from files):
   - Molecular -> drug response, `a`: rho -0.102, n 83 (output/trackA_hgcc_tmz_102/result.json; paper/trackA_hgcc_tmz_102_result_2026-11-12.md).
   - MRI -> biology, used for `b` and for the MRI-molecular link: combined rho 0.114, n 49 (output/integration_methoddev/confirmation.json).
   - SE for each = 1/sqrt(n-3) (Fisher z scale).
2. EARF ambiguity set (primary): `credible_box_model_set` � a, b and the MRI-molecular link rho each at (low, estimate, high) of its 90% Fisher-z interval; d solved from rho (clipped to +-0.95); times the 5 q levels = 135 models, same size as the current grid. Values: a in {-0.279, -0.102, 0.081}, b in {-0.127, 0.114, 0.343}.
   Diagnostic only (not used by the controller): on the coarse 135-model grid, the 90%-weight subset below keeps only a=b=d=0, so the coarse grid cannot express the evidence; this is why the credible box is used.
2b. Evidence weight of each of the 135 grid models m=(q,a,b,d):
   w(m) = N(atanh(a); atanh(a_hat), se_a) * N(atanh(b); atanh(b_hat), se_b) * N(atanh(rho_m); atanh(rho_hat), se_b), where rho_m = a*b + sqrt(1-a^2)sqrt(1-b^2)d (standard_normal_marginal_signature). Coefficients 0 and +-0.5 map through atanh unchanged in rank. q is uniform.
3. EARF = `SignalAmbiguityBeliefMPC` over the credible-box set (both channels). No other setting changes.
4. Code: `src/integration_earf.py` (`credible_box_model_set`, diagnostic `evidence_weights`/`evidence_anchored_model_set`, world sampler `sample_evidence_world`); arm `earf_minimax_mpc` in `src/211_trackc_policy_benchmark.py` when `earf_model_parameters` is passed (default None, so the Track C primary run is unchanged).

## Test (one locked test, run once on the Spark run)
- World family: `evidence_anchored` — each world draws true (a, b) from the Fisher-z normal anchors above, d such that rho matches a draw from the MRI-link anchor (clipped to [-0.95, 0.95]); q from the primary sampling levels.
- Contrasts (same patients, seeds, prior, controller code):
  - `earf_vs_full_minimax`: earf_minimax_mpc minus finite_set_minimax_mpc
  - `earf_vs_no_signal`: earf_minimax_mpc minus no_signal_mpc
- Rule: the frozen Track C rule, unchanged (GE-betting, delta 0.025 per family, both family lower bounds > 0, half-widths <= 14 days). EARF positive iff both contrasts positive.
- Secondary, descriptive: the same two contrasts at the unchanged primary anchor (0.5, 0.5, 0.5), to show the cost when signals are truly strong.

## Weak points (state with every result)
- Simulation only; no patient benefit. Sign conventions mapping rho to loadings are assumptions.
- Partly circular: worlds and the ambiguity set use the same anchors. The primary-anchor secondary analysis is the check against this.
- The anchors are proxies: MRI->MES is not MRI->drug resistance; HGCC is cell lines.
- Data-driven ambiguity sets are established in robust optimization (e.g. Delage & Ye 2010). The new part is the application: anchoring a GBM treatment controller's signal ambiguity on failed real-data checks, tested with a fair added-value design. Full prior-art search not done.
- Track C precision gate (14 days) is unchanged and was infeasible at current seed capacity; this test inherits that gate.
