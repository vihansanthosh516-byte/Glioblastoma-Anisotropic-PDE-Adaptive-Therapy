# Track C paired-signal action-sensitivity protocol (2026-11-06)

**Status: frozen Stage 2 development diagnostic.** This checks whether generated paired signals change the finite-set minimax controller's first plan in a small set of synthetic states. It does not score treatment outcomes, estimate patient benefit, or test clinical value. It changes no patient-data scope and makes no method-selection rule.

**Pre-run record correction:** An initial manifest was written before the final source-hash and metadata checks were added. No controller plan was computed from it. The initial manifest remains at `output/trackc_signal_action_sensitivity_20261106/manifest.json` and must not be used. The active manifest and result paths are under `output/trackc_signal_action_sensitivity_20261106_v2/`; commit the final manifest and seed allocation before running.

## Question

Across four synthetic virtual-patient states, do repeated paired molecular-like and MRI-like signal draws change the three-cycle action pattern selected by the full 81-model minimax controller? How often do its plans match the model set that fixes residual dependence `d=0` while retaining the same `q/a/b` support?

## Frozen design

- **Synthetic states:** one patient from each of four independent scenario seeds: `20262000`, `20262001`, `20262002`, and `20262003`.
- **Signal draws:** 16 independent paired-signal draws for each state, using seeds `20262004–20262067`. The assumed generator parameters will be `q=0.7`, `a=0.75`, `b=0.75`, and `d=0.75`. These are simulation assumptions, not estimates from Track A or Track B. `d=0.75` lies outside the controller's tested `[-0.5, 0.5]` grid.
- **Policy inputs:** each controller receives only the generated binary molecular-like report, continuous MRI-like score, and exact `log(B0)` baseline reading. It receives no latent `f0`, resistance score, signal truth, or later outcomes.
- **Controllers:** the full set has 81 combinations with `q ∈ {0.5, 0.7, 0.9}` and `a,b,d ∈ {-0.5, 0, 0.5}`. The paired control uses the same `q/a/b` values and fixes `d=0`, for 27 models.
- **Planning settings:** 128 posterior particles, 64 prediction particles, `scan_every=3`, and a three-cycle forecast. Filter initialization seed `930001001` and prediction seed `930001002` stay fixed across all signal replicates and both model sets. This pairs random streams so changes reflect signal inputs and model-set membership.
- **Measured outputs:** first selected three-cycle plan, worst-case planning-value margin, elapsed planning time, signal values, and the paired full-versus-`d=0` plan agreement. No treatment trajectory or outcome will be simulated.

## Interpretation

Report each model set's plan counts by synthetic state and exact paired agreement. No pass threshold will be applied, and no setting will be selected. Plan changes show only decision sensitivity in these test cases. Unchanged plans suggest limited action sensitivity under these settings. Neither result establishes that the signals are biologically valid or that a treatment plan helps.

## Seeds and reproducibility

- Scenario/signal seeds: `20262000–20262067`, checked against the Stage 2 development registry before preparation.
- Policy seeds: `930001001–930001002`, checked against the registry before preparation.
- No Stage 3 sealed seeds, patient MRI, masks, or real treatment data will be used.
- The pre-run manifest will hash this protocol and the diagnostic/controller/simulator sources. The result will preserve all 128 rows and summarize them without hiding ties, stable actions, or disagreement.

## Limits and next gate

The diagnostic uses four synthetic states, one assumed signal-generating model, and exact baseline readings without burden measurement noise. It is a sensitivity check, not an outcome comparison. A multi-world policy-outcome benchmark remains necessary before judging whether the method improves simulated utility. The signal-to-resistance relationship remains an assumption because Track A, Track B, and Track C patient identities are not linked.
