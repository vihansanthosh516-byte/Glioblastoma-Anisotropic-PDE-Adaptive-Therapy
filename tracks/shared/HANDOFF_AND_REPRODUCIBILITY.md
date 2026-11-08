# Independent handoff and reproducibility guide

## What a new researcher can and cannot reproduce now

The scientific question, simulation equations, signal-generating model, controller contrast, progression outcomes, and reporting rule are described in the master plan and [`configs/trackABC_integrated_v2.json`](../../configs/trackABC_integrated_v2.json). The project map identifies the existing implementation, tests, and records for each track.

The integrated study is **not yet an executable or frozen evaluation package**. The current project has prototypes, but no research-scale driver that executes the full multi-world design. The virtual-patient cohort size, development and evaluation world counts, final materialized seed list, runtime limit, input manifest, and clean-install check still must be completed before evaluation. These are explicit gates, not numbers to guess or fill after seeing results. A new researcher can reconstruct the method and implement the planned workflow now; they cannot yet reproduce a final integrated result exactly.

## Clean-checkout order

1. Read [`START_HERE.md`](../../START_HERE.md), this guide, the integrated master plan, and the track hubs.
2. Create a clean Windows environment with Python 3.14.4. Install the pinned packages from [`requirements-win-py314.lock.txt`](requirements-win-py314.lock.txt). The lock is a snapshot of the existing project `.venv`; `uv pip check` found all 152 packages compatible in that environment. A clean install on a second machine has not yet been verified.
3. Review the exact dataset and agreement terms. Start from [`trackABC_input_manifest_v2.csv`](../../configs/trackABC_input_manifest_v2.csv), which records the A5 aggregate's SHA-256 and leaves Duke out while rights are unresolved. Add rows from [`trackABC_input_manifest_template.csv`](../../configs/trackABC_input_manifest_template.csv) only for inputs whose use is documented and allowed. Keep restricted agreements, paths, identifiers, and disallowed derived outputs outside Git.
4. Implement the research-scale driver using the settings file. Do not use archived development outcomes as evaluation data.
5. Resolve the pre-evaluation gates below, materialize and reserve the seed table, freeze the config and source hashes, and only then generate evaluation outcomes.
6. Run the frozen analysis once. Preserve all failures, missing arms, negative results, and any precision-gate failure.

## Signal generation and exact joint likelihood

For each virtual patient, map the simulated resistant fraction to a latent standard-normal score:

`zR = Φ⁻¹(clip((ln(f0) − ln(L)) / (ln(U) − ln(L)), 1e−8, 1 − 1e−8))`, with `L = 0.001 exp(−0.7)` and `U = 0.9`.

Generate independent `EB, EI ~ Normal(0,1)` and set:

`zB = a zR + sqrt(1 − a²) EB`

`zI = b zR + sqrt(1 − b²) (d EB + sqrt(1 − d²) EI)`.

The hidden molecular-like state is `C = 1[zB ≥ 0]`. The observed binary report flips `C` with probability `1 − q`. The MRI-like baseline score is `zI`, or a fixed affine transformation from an allowed history-only template. If that template is not permitted or unavailable, use the standard-normal score. These are simulated signals; `a`, `b`, `d`, and `q` are assumptions, not estimates of patient-level biology.

For a particle with score `zR`, and an observed MRI-like score `zI`, compute `v = (zI − b zR) / sqrt(1 − b²)`. Then:

`pB = Φ((a zR + sqrt(1 − a²) d v) / sqrt((1 − a²)(1 − d²)))`

`pReportPositive = q pB + (1 − q)(1 − pB)`.

The joint signal likelihood is the Gaussian density `Normal(zI; b zR, 1 − b²)` multiplied by `pReportPositive` for a positive report and by `1 − pReportPositive` for a negative report. If `q = 0.5`, the molecular report contributes no information. The implementation is in `src/trackc_integrated_signal_model.py` and must be checked against this equation before it is used in the integrated driver.

## Particle updates and policy loop

1. Draw the declared parameter particles from the fixed prior in the master plan, set `S0 = (1 − f0)B0`, `R0 = f0 B0`, and assign equal particle weights.
2. For each candidate signal model, clone the same initial prior particles and add its joint molecular/MRI log likelihood. Use 500 particles per model.
3. At each visible scan, for each particle with burden `Bᵢ = Sᵢ + Rᵢ`, add `log Normal(y; ln(max(Bᵢ, 1e−300)), σ²)` to its log weight. Normalize with a stable log-sum-exp calculation. A missing scan does not update weights. If effective sample size `1 / Σwᵢ²` is below half the particle count, use the systematic resample-move step in `src/trackc_sim.py` (`h = 0.90`): systematic resampling, covariance `V = weighted_cov(z) + 1e−10 I`, proposal `z' = h z_resampled + (1−h) z_bar + chol((1−h²)V) ε`, then deterministic state re-propagation under known actions.
4. For each of the eight three-cycle treatment/no-treatment patterns, draw 16 common stratified posterior particles per model. Repeat the pattern twice to forecast six cycles (168 days) using common observation-noise and missing-scan draws. Score mean utility under every model and choose the pattern with the largest minimum model-specific mean utility. Resolve ties by the fixed order that favors fewer treatment courses.
5. Apply the first three-cycle block, then replan at the next scheduled scan. Each treatment course applies the modeled treatment for five days; no course exceeds the total budget of 12.
6. Record both latent progression times, the observed-scan treatment-stop history, course count, censoring, model-level values, and all seed keys for every virtual patient and world.

The exact ODE equations, priors, scan schedule, stop rule, utility, primary contrasts, bootstrap count, and pass rule are in the master plan and machine-readable config. The resistant fitness-cost parameter `c` will be drawn per patient from `Uniform(0, 0.25)`; it is not a second world-level axis. The existing implementation modules are prototypes; they are not yet joined into the integrated execution driver.

## World and patient sampling

Generate continuous world axes with `scipy.stats.qmc.LatinHypercube(scramble=True)` using the world's named seed. Use the dimension order in the machine-readable config, then map each unit-cube coordinate through its declared linear or log-uniform range. Assign the nine `(mu, nu)` combinations in lexicographic order, repeat them, randomly permute the repeated list with the world-parameter stream, and take the first `N` worlds; counts across combinations will differ by at most one. Generate patient parameters with a separate scrambled Latin hypercube in the configured dimension order. For mixture worlds, assign response class with its own stream and then apply the responder or nonresponder `m` range. Assign `sigma` by repeating `[0.05, 0.10, 0.20]` across sorted patient IDs, permuting once with the observation-noise stream, and taking the assigned values; group counts will differ by at most one. Every exact cohort size and LHS seed is written in the materialized pre-run seed table.

## Seed table

[`trackABC_seed_streams_v2.csv`](../../configs/trackABC_seed_streams_v2.csv) defines the names and order-independent key for each random stream. The integer seed is the little-endian value of the first eight bytes of `SHA256(UTF8("trackABC-integrated-v2|stage|family|world_id|stream_name"))`; initialize `numpy.random.Generator(numpy.random.PCG64DXSM(seed))`. Before use, materialize every numeric seed into a dated table, check duplicates and collisions against `src/trackc_seed_registry.py`, and add the appropriate reservation. The CSV currently contains role templates, not allocated or consumed seeds.

Development, selection, calibration, evaluation, molecular signal, MRI signal, observation noise, each policy, and world bootstrap must use distinct stream names. Within a world, all policies must receive the same patients, latent parameter draws, baseline signals, scan schedule, and common noise draws wherever a paired comparison requires them. Never use an existing Stage 3 sealed reservation for development.

## Input checksums and use rights

The integrated input manifest at [`configs/trackABC_input_manifest_v2.csv`](../../configs/trackABC_input_manifest_v2.csv) retains the A5 aggregate file hash for provenance only; observed A5 metrics do not set q_signal or otherwise parameterize the integrated model. It does not include a Duke-derived summary. The run manifest will hash only the exact allowed aggregate files that the run reads. Record SHA-256, file size, source record, license/permission basis, content type, and allowed use. Dataset rights and article rights are separate. The Duke agreement and record terms govern every Duke-derived quantity; if derived-summary rights are unclear, omit the Duke summary and use the standard-normal MRI-like fallback. Do not hash, copy, or reopen raw patient scans, masks, ADC maps, or private patient records for this phase.

## Items to lock before evaluation

- Select virtual-patient cohort sizes and the number of development worlds using a written outcome-blind timing and precision rule; preserve the planning data separately from evaluation data.
- State the formula and thresholds for that planning rule, including the maximum evaluation-world count and the action if the 14-day interval half-width target cannot be met.
- State the exact Latin-hypercube algorithm and the deterministic assignment order for balanced observation-noise levels.
- Materialize and reserve the numeric seed list, complete the integrated driver, and record hashes for source, config, input summaries, and environment.
- Verify that the Duke-derived aggregate, if any, may be used under the written agreement; otherwise keep the fallback.
- Set and record a runtime cap from timing-only measurements. If the gate fails, stop before evaluation or report feasibility only.

## Do not interpret as

The integrated study is a synthetic comparison under declared assumptions. It cannot establish a clinical treatment benefit, a resistance biomarker, patient-level fusion, or superiority over all published methods. Track A and Track B findings remain prior exploratory evidence; Track C development results remain exposed development evidence.
