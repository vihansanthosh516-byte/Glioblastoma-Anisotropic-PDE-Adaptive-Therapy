# Integrated signal-to-resistance bridge (proposed draft, 2026-11-05)

## Status

This is a proposed synthetic-method amendment. It is not adopted, frozen, or authorization to run a held-out evaluation. It changes no patient-data scope. Unit tests will check implementation behavior only; they will not establish a positive treatment result or biological validity.

## Problem being fixed

The earlier signal generator drew its resistance score independently of the simulator's initial resistant fraction `f0`. A policy could not use those signals to update its belief about the state that drives simulated resistance. This draft connects the signal generator to `f0` through one fixed, explicit transform.

## Fixed resistance score

For every virtual patient, let `f0` be the initial resistant-cell fraction drawn by the existing calibration family. Define the global reference bounds as:

`f_min = 0.001 * exp(-0.7)` and `f_max = 0.9`.

These bounds cover the union of the assumed single and mixture calibration-family ranges. They are global across worlds; they are not refit within a hidden world. Define:

`u = clip(log(f0 / f_min) / log(f_max / f_min), 1e-8, 1 - 1e-8)`

`Z_R = Phi_inverse(u)`.

Thus the same `f0` always maps to the same resistance score. This is a mathematical linking assumption, not a measured molecular-to-resistance relationship.

## Signal-generating model

For independent standard-normal errors `E_B` and `E_I`, define:

`Z_B = a * Z_R + sqrt(1 - a^2) * E_B`

`Z_I = b * Z_R + sqrt(1 - b^2) * (d * E_B + sqrt(1 - d^2) * E_I)`

The synthetic molecular state is `C = 1[Z_B >= 0]`. A hypothetical binary report `M` reports `C` correctly with probability `q`; this uses symmetric sensitivity and specificity. `Z_I` is the synthetic baseline MRI/ADC score. If a Track B score template is supplied later, its empirical marginal will be fitted on permitted training patients only and applied by a monotone rank transform.

`a`, `b`, `d`, and `q` are analyst-chosen scenario assumptions. A5 did not measure resistance, did not establish this mapping, and did not validate a clinical assay. No patient-level Track A or Duke rows will be joined.

## Policy likelihood and comparison

The signal-informed belief MPC will receive only `M`, the MRI/ADC score, and one fixed signal-likelihood setting selected before evaluation. It will not receive `f0`, `Z_R`, `Z_B`, or any other simulated truth. The fixed-prior primary likelihood will use the predeclared moderate setting; sensitivity worlds will vary their data-generating `a`, `b`, `d`, and `q` while the policy likelihood stays fixed. This avoids giving the policy each evaluation world's true signal mechanism.

At initialization, the particle filter will add the joint log likelihood `log p(M, Z_I | f0)` to each particle's log weight. Later longitudinal burden scans will update the same particle filter under the existing scan likelihood. The no-signal, molecular-only, MRI/ADC-only, and two-signal policies will use the same broad `f0` prior, virtual patients, random streams, treatment limits, and evaluation code. A no-signal policy will apply no initial likelihood. Signal-shuffle controls will preserve signal marginals while breaking patient alignment.

If an empirical MRI template is used, each held-out score will be mapped through a mid-rank Gaussian transform fitted only on permitted training data. The resulting update is a rank-Gaussian pseudo-likelihood for relative particle weights; it is not a calibrated probability model for the raw MRI score. The template and transform will be fixed before held-out evaluation, and their limitations will be reported.

## Synthetic implementation checks

Before any development outcome simulation, tests will check:

1. `f0` maps monotonically through one global reference transform and covers every assumed calibration-family bound.
2. Generated signals follow the declared sign of `a` and `b`; the joint residual dependence follows `d`.
3. At `a = b = 0`, signal likelihood is constant over `f0` particles, even when `d` is nonzero.
4. Reversing `a` or `b` reverses the corresponding likelihood ordering.
5. The empirical MRI rank transform gives finite likelihoods and uses only its supplied template.
6. The policy's initialized weights and actions do not change when hidden true parameters change while supplied observations, visible scan precision, seed, and policy assumptions stay fixed.
7. The no-signal configuration matches the ordinary belief MPC under the same broad prior.
8. Invalid shapes, non-finite observations, and invalid parameter values fail clearly.

Passing these checks will establish only code/model consistency. It will not establish that the assumptions describe glioblastoma or that fusion improves therapy outcomes.

## Remaining gates

- The marker benchmark selector now implements the integrated plan's patient-wise utility and exact tie rule, with focused unit tests. It is not yet wired into the shared development benchmark; that integration remains open.
- Freeze the primary signal likelihood, reliability values, signal-shuffle protocol, seeds, and development/held-out split before outcome simulations.
- Decide how approved Track B-derived MRI/ADC summaries enter the synthetic marginal; no Duke voxel data are read by this draft or its unit tests.
- Keep robust-policy and Q-learning comparisons separate until each method's inputs and training information match this protocol.
- Update the master and short research plans only after this amendment is reviewed and adopted.
