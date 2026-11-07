# Track C minimax action-sampling stability protocol (2026-11-06)

**Status: frozen development diagnostic; not yet run.** This test checks Monte Carlo stability of controller choices. It does not test treatment benefit, patient benefit, or biological validity. It does not use patient data or Stage 3 seeds.

## Question

Do the finite-set minimax controller's three-cycle treatment plans remain unchanged as the posterior prediction sample count increases from 16 to 32 to 64, under the same posterior filters and synthetic baseline state?

## Frozen design

- **Controller:** `SignalAmbiguityBeliefMPC`, with a three-cycle forecast and a three-cycle scan interval. Posterior filter size will be 128 particles. The simulator's default course budget and five-day course penalty will be used.
- **Synthetic state:** one virtual patient from scenario seed `20261990`. The baseline burden reading will be exact `log(B0)`. No reading noise will be added.
- **Signal anchors:** four predeclared input pairs will be tested: `(molecular=True, MRI=+0.75)`, `(True, -0.75)`, `(False, +0.75)`, and `(False, -0.75)`. These are analyst-chosen test inputs, not patient measurements.
- **Model sets:** the full 81-point grid uses `q ∈ {0.5, 0.7, 0.9}` and `a,b,d ∈ {-0.5, 0, 0.5}`. The matched comparator uses the same `q/a/b` values and fixes `d=0`, for 27 models.
- **Posterior filters:** every particle-count condition for a given anchor and model set will use the same 128-particle filter initialization, policy seed `930000004`, patient, and baseline burden reading.
- **Prediction draws:** each condition will use particle counts `{16, 32, 64}` and prediction seeds `{930000001, 930000002, 930000003}`. Keyed stratified draws will be shared across signal models within each condition.
- **Measured output:** the selected three-cycle treatment pattern, each pattern's minimum model-specific predicted utility, the margin between the best and second-best patterns, and elapsed time. No TTP, survival, course outcome, or treatment-benefit metric will be computed.

## Interpretation rule

The primary diagnostic is exact plan agreement between adjacent prediction-particle counts, paired by signal anchor, model set, and prediction seed. Any plan change from 32 to 64 will mark that setting unstable for the tested anchors. Full-grid versus matched-`d=0` plan agreement will be reported separately. A stable plan in these four anchors will support only numerical stability in these test cases; it will not establish general stability, useful treatment effects, or a valid signal-to-resistance relationship. No model set or setting will be selected from treatment outcomes.

## Seed allocation

The scenario seed will be `20261990` from the Stage 2 development range. Policy RNG seeds will be `930000001` through `930000004`. The seed registry will record the run after the pre-run manifest is written. No noise-key stream is used.

## Reproducibility and limits

The pre-run manifest will hash the diagnostic code, controller code, simulator code, signal-model code, and this protocol. The result will contain action plans and predicted planning values only. One synthetic patient state, no observation noise, and analyst-chosen signal anchors make this a narrow numerical diagnostic. The simulator's signal-to-resistance mapping remains an assumption.
