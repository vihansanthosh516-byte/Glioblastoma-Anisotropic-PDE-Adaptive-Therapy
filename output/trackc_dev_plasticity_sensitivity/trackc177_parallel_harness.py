import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import trackc_sim as tc

MU_NU_RATES = (0.0, 0.01, 0.04)
SCENARIO_SEED = 20262200
NOISE_SEED = 1_620_000
POLICIES = {
    "front_loaded_q1": ("periodic", 1, 0),
    "raw_threshold_0.6": ("raw", 0.6, 0.6),
    "wide_containment": ("riskf", 0.1, 2, 0),
    "holiday_max_run_3": ("riskh", 0.1, 2, 0, 3),
}


def noise_stream_digest(seed, n, sigma):
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal((n, tc.N_CYCLES)) * sigma[:, None]
    return hashlib.sha256(np.ascontiguousarray(noise).tobytes()).hexdigest()


def paired_bootstrap_mean(diff, n_boot, seed):
    n_rep, n_scen = diff.shape
    point = float(diff.mean())
    rng = np.random.default_rng(seed)
    estimates = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        s_idx = rng.integers(0, n_scen, n_scen)
        r_idx = rng.integers(0, n_rep, (n_rep, n_scen))
        estimates[b] = diff[r_idx, np.broadcast_to(s_idx, (n_rep, n_scen))].mean()
    return point, np.percentile(estimates, [2.5, 97.5]).tolist()


def compute_cell(args):
    mu, nu, n, replicates, n_boot, particles = args
    scen = {k: np.asarray(v).copy() for k, v in tc.draw_scenarios(n, SCENARIO_SEED).items()}
    scen["mu"] = np.full(n, mu, dtype=float)
    scen["nu"] = np.full(n, nu, dtype=float)
    seeds = [NOISE_SEED + r for r in range(replicates)]
    outcomes = {}
    cell_results = {}
    for name, spec in POLICIES.items():
        policy = tc.make_policy(spec)
        if hasattr(policy, "n_particles"):
            policy.n_particles = particles
        runs = [tc.simulate(scen, policy, seed) for seed in seeds]
        ttp_reps = np.stack([x["ttp"] for x in runs])
        course_reps = np.stack([x["courses"] for x in runs])
        outcomes[name] = {
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
    front = np.asarray(outcomes["front_loaded_q1"]["ttp_days_replicate_by_scenario"])
    for name in POLICIES:
        if name == "front_loaded_q1":
            continue
        candidate = np.asarray(outcomes[name]["ttp_days_replicate_by_scenario"])
        mean_diff, ci = paired_bootstrap_mean(candidate - front, n_boot, seed=177)
        cell_results[name]["vs_front_loaded_mean_diff_days"] = mean_diff
        cell_results[name]["vs_front_loaded_mean_diff_ci95_days"] = ci
        cell_results[name]["vs_front_loaded_cvar10_scenario_mean_diff_days"] = float(
            tc.cvar_scenario_mean(candidate) - tc.cvar_scenario_mean(front)
        )
    return f"mu={mu:g},nu={nu:g}", cell_results, outcomes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--replicates", type=int, required=True)
    ap.add_argument("--bootstrap-draws", type=int, required=True)
    ap.add_argument("--particles", type=int, required=True)
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    start = time.time()
    jobs = [(mu, nu, a.n, a.replicates, a.bootstrap_draws, a.particles)
            for mu in MU_NU_RATES for nu in MU_NU_RATES]
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        rows = list(ex.map(compute_cell, jobs))
    cells, outcomes = {}, {}
    for key, summary, scenario_data in rows:
        cells[key] = summary
        outcomes[key] = scenario_data
    base = tc.draw_scenarios(a.n, SCENARIO_SEED)
    seeds = [NOISE_SEED + r for r in range(a.replicates)]
    result = {
        "script": "177_trackc_dev_plasticity_sensitivity",
        "label": "SIMULATION ONLY; development-only; exploratory; assumed mu/nu stress values",
        "protocol_status": "new development screen; no sealed hypotheses or test seeds; no existing protocol changed",
        "n_scenarios": a.n,
        "noise_replicates": a.replicates,
        "bootstrap_draws": a.bootstrap_draws,
        "filter_particles": a.particles,
        "scenario_seed": SCENARIO_SEED,
        "noise_seeds": seeds,
        "rates_per_day_assumed": list(MU_NU_RATES),
        "base_scenarios_shared_across_all_cells": True,
        "observation_noise": {
            "shared_exactly_across_policies_and_rate_cells": True,
            "reason": "tc.simulate draws the full observation-noise matrix first from the replicate seed and same sigma array; rate cells do not alter sigma or draw shape",
            "sha256_by_seed": {str(seed): noise_stream_digest(seed, a.n, base["sigma"]) for seed in seeds},
            "limitation": "filter random draws start after observation noise and share initial RNG stream, but posterior-dependent resampling can make later filter states diverge; this is not common random numbers for every internal filter decision",
        },
        "policies": {k: list(v) for k, v in POLICIES.items()},
        "results": cells,
        "per_scenario_outcomes": outcomes,
        "runtime_s": time.time() - start,
    }
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "dev_sensitivity_results.json").write_text(json.dumps(result, indent=2))
    print(f"wrote {a.out / 'dev_sensitivity_results.json'}")
    print(f"cells={len(cells)}; n={a.n}; reps={a.replicates}; runtime={result['runtime_s']:.1f}s")


if __name__ == "__main__":
    main()
