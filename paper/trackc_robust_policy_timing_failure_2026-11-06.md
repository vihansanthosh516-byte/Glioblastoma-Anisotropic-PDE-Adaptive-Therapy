# Robust policy timing pilot: stopped at input guard

**Status:** failed before policy simulation; no policy outcomes saved or inspected.

The frozen two-family timing pilot stopped during its first (`single`) family. The runner rejected the generated evaluation scenarios because their `f0` values fell outside `PolicyBenchmarkSettings.signal_f0_bounds` (`0.001`–`0.1`). The scenario generator shifts the `f0` range by calibration world, so its generated values can exceed those default bounds. The exception was `ValueError: evaluation_scenarios contains f0 values outside the declared signal_f0_bounds`.

The run reached the input guard before any policy simulation. It did not run the `mixture` family. It saved no policy score, utility, endpoint, course count, or ranking. No patient data, MRI, mask, or restricted file was read. The exact environment, source hashes, protocol hash, seed allocation, and failure are preserved in `output/trackc_robust_policy_timing_pilot_20261106/manifest.json`.

This attempt used the `single` world and input-generation seeds 20261930–20261935. It consumed no policy RNG seed. These scenario seeds are registered as development-used to prevent reuse. The original protocol's mixture seeds and unconsumed policy seeds remain unused. Any follow-up timing pilot must have its own dated amendment, use fresh seeds, and set the signal-model bounds to a fixed envelope that covers the calibration generator's declared `f0` ranges. The failed attempt is not a performance result.
