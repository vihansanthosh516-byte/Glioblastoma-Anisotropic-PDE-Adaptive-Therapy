"""Finite-sample recovery check for the assumed Track C paired-signal model.

This development diagnostic uses only synthetic signal observations. It estimates
no patient parameters and evaluates no treatment policy. The protocol and fixed
settings are in paper/trackc_finite_sample_recovery_protocol_2026-11-06.md.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from trackc_integrated_signal_model import (
    SignalModelParameters,
    SignalObservations,
    generate_integrated_signals,
    integrated_signal_log_likelihood,
)


ROOT = Path(__file__).resolve().parent.parent


def _load_sibling(name: str, filename: str):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {filename}: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_GRID = _load_sibling("trackc_signal_fisher_grid_for_recovery", "212_trackc_signal_fisher_grid.py")
_DEPENDENCE = _load_sibling("trackc_dependence_identifiability_for_recovery", "210_trackc_dependence_identifiability.py")
_SEED_REGISTRY = _load_sibling("trackc_seed_registry_for_recovery", "trackc_seed_registry.py")


F0_RANGE_NAMES = ("calibration_low_shift", "default", "calibration_high_shift")
F0_RANGES = {name: _GRID.F0_RANGES[name] for name in F0_RANGE_NAMES}
TRUTH_ANCHORS = (
    (0.7, 0.5, 0.5, -0.5),
    (0.7, 0.5, 0.5, 0.0),
    (0.7, 0.5, 0.5, 0.5),
    (0.9, 0.5, 0.5, 0.5),
    (0.5, 0.5, 0.5, 0.5),
    (0.7, 0.0, 0.5, 0.5),
    (0.7, 0.5, 0.0, 0.5),
)
SAMPLE_SIZES = (18, 28, 100, 500)
ROOT_SEEDS = tuple(range(20261880, 20261905))
N_F0_NODES = 256
PROFILE_LOG_LIKELIHOOD_DROP = 2.0
TIE_TOLERANCE = 1e-8


def parameter_grid():
    """Return the fixed 81 candidates from the Fisher stress-test protocol."""
    return _GRID.parameter_grid()


def _quadrature(bounds, n_nodes=N_F0_NODES):
    """Return f0 nodes and normalized weights for the log-uniform prior."""
    return _DEPENDENCE.log_uniform_f0_quadrature(bounds, n_nodes=n_nodes)


def marginal_log_likelihoods(observations, f0_bounds, *, n_nodes=N_F0_NODES):
    """Score every candidate after integrating out each unobserved patient's f0."""
    reports = np.asarray(observations.molecular_report, dtype=bool)
    scores = np.asarray(observations.mri_score, dtype=float)
    if reports.ndim != 1 or scores.shape != reports.shape or reports.size < 1:
        raise ValueError("observations must contain equal, nonempty one-dimensional arrays")
    if not np.isfinite(scores).all():
        raise ValueError("MRI scores must be finite")

    f0_nodes, weights = _quadrature(f0_bounds, n_nodes=n_nodes)
    particles = np.broadcast_to(f0_nodes, (reports.size, f0_nodes.size))
    log_weights = np.log(weights)[None, :]
    values = []
    for candidate in parameter_grid():
        conditional = integrated_signal_log_likelihood(
            particles,
            SignalObservations(reports, scores),
            SignalModelParameters(*candidate),
        )
        marginal = logsumexp(conditional + log_weights, axis=1)
        value = float(np.sum(marginal))
        if not np.isfinite(value):
            raise ArithmeticError(f"non-finite marginal log likelihood at {candidate}")
        values.append(value)
    return np.asarray(values, dtype=float)


def _seeded_dataset(seed_root, range_index, truth_index, size_index, truth, bounds, n):
    """Generate observations with keyed streams; keep latent f0 out of the estimator."""
    seed_sequence = np.random.SeedSequence(
        [int(seed_root), int(range_index), int(truth_index), int(size_index)]
    )
    f0_sequence, signal_sequence = seed_sequence.spawn(2)
    f0_rng = np.random.default_rng(f0_sequence)
    signal_seed = int(signal_sequence.generate_state(1, dtype=np.uint32)[0])
    low, high = map(float, bounds)
    f0 = np.exp(f0_rng.uniform(np.log(low), np.log(high), size=int(n)))
    observations, _truth = generate_integrated_signals(
        int(n), q=truth[0], a=truth[1], b=truth[2], d=truth[3],
        seed=signal_seed, f0=f0,
    )
    return observations


def evaluate_observations(observations, truth, f0_bounds, *, n_nodes=N_F0_NODES):
    """Return transparent grid-MLE and descriptive near-maximum profile diagnostics."""
    candidates = parameter_grid()
    likelihoods = marginal_log_likelihoods(observations, f0_bounds, n_nodes=n_nodes)
    if likelihoods.shape != (len(candidates),):
        raise ArithmeticError("likelihood vector does not match candidate grid")

    maximum = float(np.max(likelihoods))
    mle_indices = np.flatnonzero(maximum - likelihoods <= TIE_TOLERANCE)
    profile_indices = np.flatnonzero(maximum - likelihoods <= PROFILE_LOG_LIKELIHOOD_DROP)
    mle_points = [candidates[int(i)] for i in mle_indices]
    profile_points = [candidates[int(i)] for i in profile_indices]
    d_hat = float(np.mean([point[3] for point in mle_points]))
    truth_tuple = tuple(map(float, truth))
    truth_index = candidates.index(truth_tuple)
    truth_log_likelihood = float(likelihoods[truth_index])
    truth_in_mle_ties = bool(truth_index in mle_indices)
    d_profile = [point[3] for point in profile_points]
    # When q=0.5, the molecular report is a coin flip and d is structurally
    # unidentifiable; do not score its sign as a recovery target.
    sign_identifiable = truth_tuple[3] != 0.0 and truth_tuple[0] != 0.5
    return {
        "mle_parameters": list(map(float, mle_points[0])),
        "mle_tie_count": int(len(mle_points)),
        "truth_in_mle_ties": truth_in_mle_ties,
        "truth_log_likelihood_minus_max": truth_log_likelihood - maximum,
        "d_hat_tie_mean": d_hat,
        "d_absolute_error": abs(d_hat - truth_tuple[3]),
        "d_sign_correct": (bool(np.sign(d_hat) == np.sign(truth_tuple[3])) if sign_identifiable else None),
        "d_profile_min_delta_2": float(min(d_profile)),
        "d_profile_max_delta_2": float(max(d_profile)),
        "d_profile_width_delta_2": float(max(d_profile) - min(d_profile)),
        "truth_d_in_profile_delta_2": bool(
            min(d_profile) - TIE_TOLERANCE <= truth_tuple[3] <= max(d_profile) + TIE_TOLERANCE
        ),
        "profile_candidate_count_delta_2": int(len(profile_points)),
    }


def _summarize(rows):
    out = {}
    keys = sorted({(r["f0_range"], r["truth_label"], r["n"]) for r in rows})
    for f0_name, truth_label, n in keys:
        group = [r for r in rows if (r["f0_range"], r["truth_label"], r["n"]) == (f0_name, truth_label, n)]
        key = f"{f0_name}|{truth_label}|n={n}"
        d_errors = np.asarray([r["d_absolute_error"] for r in group], dtype=float)
        d_widths = np.asarray([r["d_profile_width_delta_2"] for r in group], dtype=float)
        truth_tie = np.asarray([r["truth_in_mle_ties"] for r in group], dtype=float)
        truth_profile = np.asarray([r["truth_d_in_profile_delta_2"] for r in group], dtype=float)
        sign_rows = [
            r for r in group
            if r["d_sign_correct"] is not None and r["truth_parameters"][0] != 0.5
        ]
        out[key] = {
            "n_replicates": len(group),
            "truth_in_mle_ties_rate": float(np.mean(truth_tie)),
            "median_absolute_d_error": float(np.median(d_errors)),
            "d_sign_recovery_rate_nonzero_truths": (
                float(np.mean([r["d_sign_correct"] for r in sign_rows])) if sign_rows else None
            ),
            "median_profile_d_width_delta_2": float(np.median(d_widths)),
            "truth_d_in_profile_delta_2_rate": float(np.mean(truth_profile)),
        }
    return out


def build_diagnostic(
    *,
    n_replicates=25,
    sample_sizes=SAMPLE_SIZES,
    n_nodes=N_F0_NODES,
    seed_roots=ROOT_SEEDS,
    check_registered_seeds=True,
):
    """Run the fixed development design; arguments may be reduced only by smoke tests."""
    seed_roots = tuple(int(seed) for seed in seed_roots)
    if isinstance(n_replicates, bool) or int(n_replicates) != n_replicates or not 1 <= n_replicates <= len(seed_roots):
        raise ValueError("n_replicates must be an integer from 1 to the 25 registered roots")
    if len(set(seed_roots)) != len(seed_roots) or any(seed < 0 for seed in seed_roots):
        raise ValueError("seed_roots must be unique non-negative integers")
    sample_sizes = tuple(int(n) for n in sample_sizes)
    if not sample_sizes or any(n < 1 for n in sample_sizes):
        raise ValueError("sample_sizes must contain positive integers")
    if isinstance(n_nodes, bool) or int(n_nodes) != n_nodes or n_nodes < 2:
        raise ValueError("n_nodes must be an integer of at least 2")

    seed_roots = seed_roots[:int(n_replicates)]
    if check_registered_seeds:
        problems = _SEED_REGISTRY.check_free(
            scenario=(seed_roots[0], seed_roots[-1]), purpose="stage2_dev"
        )
        if problems:
            raise RuntimeError(f"registered development seed check failed: {problems}")

    rows = []
    candidates = parameter_grid()
    for range_index, f0_name in enumerate(F0_RANGE_NAMES):
        bounds = F0_RANGES[f0_name]
        for truth_index, truth in enumerate(TRUTH_ANCHORS):
            truth_label = f"q{truth[0]}_a{truth[1]}_b{truth[2]}_d{truth[3]}"
            for size_index, n in enumerate(sample_sizes):
                for seed_root in seed_roots:
                    observations = _seeded_dataset(
                        seed_root, range_index, truth_index, size_index, truth, bounds, n
                    )
                    result = evaluate_observations(observations, truth, bounds, n_nodes=n_nodes)
                    rows.append({
                        "f0_range": f0_name,
                        "f0_bounds": list(map(float, bounds)),
                        "truth_label": truth_label,
                        "truth_parameters": list(map(float, truth)),
                        "n": int(n),
                        "seed_root": int(seed_root),
                        **result,
                    })

    expected_rows = 3 * len(TRUTH_ANCHORS) * len(sample_sizes) * int(n_replicates)
    if len(rows) != expected_rows:
        raise ArithmeticError(f"expected {expected_rows} rows, found {len(rows)}")
    return {
        "scope": "synthetic finite-sample signal-parameter recovery; no patient data or treatment outcomes",
        "protocol": {
            "f0_ranges": {name: list(map(float, F0_RANGES[name])) for name in F0_RANGE_NAMES},
            "truth_anchors": [list(map(float, row)) for row in TRUTH_ANCHORS],
            "candidate_grid": [list(map(float, row)) for row in candidates],
            "candidate_count": len(candidates),
            "sample_sizes": list(sample_sizes),
            "replicates_per_condition": int(n_replicates),
            "scenario_seed_roots": [int(seed_roots[0]), int(seed_roots[-1])],
            "quadrature_nodes": int(n_nodes),
            "profile_log_likelihood_drop": PROFILE_LOG_LIKELIHOOD_DROP,
            "tie_tolerance": TIE_TOLERANCE,
        },
        "summary": _summarize(rows),
        "all_replicates": rows,
    }


def main():
    print(json.dumps(build_diagnostic(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["build_diagnostic", "evaluate_observations", "marginal_log_likelihoods", "parameter_grid"]
