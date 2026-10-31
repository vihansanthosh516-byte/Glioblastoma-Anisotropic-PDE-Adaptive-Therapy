#!/usr/bin/env python3
"""Development-only paired sensitivity screen for assumed reversible tolerance.

SIMULATION ONLY. Never reads sealed outcomes or test seeds. This script does not
change any Track C protocol. It varies assumed mu (on-drug S->R) and nu
(off-drug R->S) independently, then compares fixed policies on common scenarios.

Default run is intentionally small. Use --out to choose a new output location;
the script refuses to overwrite an existing result. Reduced smoke runs must
also use a separate --out directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_sim as tc  # noqa: E402

OUT_DIR = PROJECT_ROOT / "output" / "trackc_dev_plasticity_sensitivity"
MU_NU_RATES = (0.0, 0.01, 0.04)  # per day; assumed stress values, not measured rates
SCENARIO_SEED = 20262200
NOISE_SEED = 1_620_000
POLICIES = {
    "front_loaded_q1": ("periodic", 1, 0),
    "raw_threshold_0.6": ("raw", 0.6, 0.6),
    "wide_containment": ("riskf", 0.1, 2, 0),
    "holiday_max_run_3": ("riskh", 0.1, 2, 0, 3),
}


def noise_stream_digest(seed: int, n: int, sigma: np.ndarray) -> str:
    """Fingerprint observation-noise draws made first inside tc.simulate."""
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal((n, tc.N_CYCLES)) * sigma[:, None]
    return hashlib.sha256(np.ascontiguousarray(noise).tobytes()).hexdigest()


def paired_bootstrap_mean(diff: np.ndarray, n_boot: int, seed: int) -> tuple[float, list[float]]:
    """Paired percentile interval; resample replicate indices within scenarios."""
    # diff shape: replicates x scenarios; policy pairing is already in each diff.
    n_rep, n_scen = diff.shape
    point = float(diff.mean())
    if n_scen == 0 or n_rep == 0:
        raise ValueError("empty outcome array")
    rng = np.random.default_rng(seed)
    estimates = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        s_idx = rng.integers(0, n_scen, n_scen)
        r_idx = rng.integers(0, n_rep, (n_rep, n_scen))
        estimates[b] = diff[r_idx, np.broadcast_to(s_idx, (n_rep, n_scen))].mean()
    return point, np.percentile(estimates, [2.5, 97.5]).tolist()


def run(n: int, replicates: int, n_boot: int, particles: int, out_dir: Path) -> dict:
    result_path = out_dir / "dev_sensitivity_results.json"
    if result_path.exists():
        raise FileExistsError(f"refusing to overwrite {result_path}")
    out_dir.mkdir(parents=True, exist_ok=True)

    # One underlying draw reused across every rate cell. Same replicate seeds
    # reused across policies and cells; simulate draws observation noise before
    # policy initialization, so observation-noise matrices match exactly.
    base = tc.draw_scenarios(n, SCENARIO_SEED)
    noise_seeds = [NOISE_SEED + r for r in range(replicates)]
    stream_hashes = {str(seed): noise_stream_digest(seed, n, base["sigma"])
                     for seed in noise_seeds}
    started = time.time()
    cells: dict[str, dict] = {}
    full_outcomes: dict[str, dict[str, dict[str, list]]] = {}

    for mu in MU_NU_RATES:
        for nu in MU_NU_RATES:
            cell = f"mu={mu:g},nu={nu:g}"
            scen = {k: np.asarray(v).copy() for k, v in base.items()}
            scen["mu"] = np.full(n, mu, dtype=float)
            scen["nu"] = np.full(n, nu, dtype=float)
            full_outcomes[cell] = {}
            cell_results = {}
            for name, spec in POLICIES.items():
                policy = tc.make_policy(spec)
                if hasattr(policy, "n_particles"):
                    policy.n_particles = particles
                runs = [tc.simulate(scen, policy, seed) for seed in noise_seeds]
                ttp_reps = np.stack([x["ttp"] for x in runs])
                course_reps = np.stack([x["courses"] for x in runs])
                full_outcomes[cell][name] = {
                    "ttp_days_replicate_by_scenario": ttp_reps.tolist(),
                    "course_counts_replicate_by_scenario": course_reps.tolist(),
                }
                cell_results[name] = {
                    "mean_ttp_days": float(ttp_reps.mean()),
                    "cvar10_scenario_mean_days": float(tc.cvar_scenario_mean(ttp_reps)),
                    "cvar10_joint_days": float(tc.cvar_joint(ttp_reps)),
                    "mean_courses": float(course_reps.mean()),
                    "capped_fraction": float(np.mean([x["capped"] for x in runs])),
                }

            front = np.asarray(full_outcomes[cell]["front_loaded_q1"]["ttp_days_replicate_by_scenario"])
            for name in POLICIES:
                if name == "front_loaded_q1":
                    continue
                candidate = np.asarray(full_outcomes[cell][name]["ttp_days_replicate_by_scenario"])
                diff = candidate - front
                mean_diff, ci = paired_bootstrap_mean(diff, n_boot, seed=177)
                cell_results[name]["vs_front_loaded_mean_diff_days"] = mean_diff
                cell_results[name]["vs_front_loaded_mean_diff_ci95_days"] = ci
                cell_results[name]["vs_front_loaded_cvar10_scenario_mean_diff_days"] = (
                    float(tc.cvar_scenario_mean(candidate) - tc.cvar_scenario_mean(front)))
            cells[cell] = cell_results

    result = {
        "script": "177_trackc_dev_plasticity_sensitivity",
        "label": "SIMULATION ONLY; development-only; exploratory; assumed mu/nu stress values",
        "protocol_status": "new development screen; no sealed hypotheses or test seeds; no existing protocol changed",
        "n_scenarios": n,
        "noise_replicates": replicates,
        "bootstrap_draws": n_boot,
        "filter_particles": particles,
        "scenario_seed": SCENARIO_SEED,
        "noise_seeds": noise_seeds,
        "rates_per_day_assumed": list(MU_NU_RATES),
        "base_scenarios_shared_across_all_cells": True,
        "observation_noise": {
            "shared_exactly_across_policies_and_rate_cells": True,
            "reason": "tc.simulate draws the full observation-noise matrix first from the replicate seed and same sigma array; rate cells do not alter sigma or draw shape",
            "sha256_by_seed": stream_hashes,
            "limitation": "filter random draws start after observation noise and share initial RNG stream, but posterior-dependent resampling can make later filter states diverge; this is not common random numbers for every internal filter decision",
        },
        "policies": {k: list(v) for k, v in POLICIES.items()},
        "results": cells,
        "per_scenario_outcomes": full_outcomes,
        "runtime_s": time.time() - started,
    }
    tmp_path = out_dir / f".dev_sensitivity_results.{os.getpid()}.tmp"
    with tmp_path.open("xb") as handle:
        handle.write(json.dumps(result, indent=2).encode("utf-8"))
        handle.flush()
        os.fsync(handle.fileno())
    try:
        # Same-volume hard-link creation is atomic and refuses an existing target.
        os.link(tmp_path, result_path)
    finally:
        tmp_path.unlink(missing_ok=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=36)
    parser.add_argument("--replicates", type=int, default=2)
    parser.add_argument("--bootstrap-draws", type=int, default=1000)
    parser.add_argument("--particles", type=int, default=500)
    parser.add_argument("--out", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    if args.n < 2 or args.replicates < 1 or args.bootstrap_draws < 1 or args.particles < 1:
        parser.error("n must be >= 2; replicates, bootstrap-draws, and particles must be >= 1")
    if args.out.resolve() == OUT_DIR.resolve() and (
        args.n != 36 or args.replicates != 2 or args.bootstrap_draws != 1000 or args.particles != 500
    ):
        parser.error("non-default smoke runs must use a separate --out directory")
    result = run(args.n, args.replicates, args.bootstrap_draws, args.particles, args.out)
    print(f"wrote {args.out / 'dev_sensitivity_results.json'}")
    print(f"cells={len(result['results'])}; n={args.n}; reps={args.replicates}; runtime={result['runtime_s']:.1f}s")


if __name__ == "__main__":
    main()
