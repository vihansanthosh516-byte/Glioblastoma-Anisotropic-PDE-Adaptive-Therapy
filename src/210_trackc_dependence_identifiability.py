"""Marginal paired-signal likelihood diagnostic for Track C assumptions.

This module integrates over a declared log-uniform f0 range. It can show whether
different (a, b, d) tuples induce distinguishable paired-signal distributions under a
specified f0 prior. It does not estimate patient parameters or score treatment policies.
"""
from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

from trackc_integrated_signal_model import (
    F0_REFERENCE_BOUNDS,
    SignalModelParameters,
    SignalObservations,
    integrated_signal_log_likelihood,
    standard_normal_marginal_signature,
)


def log_uniform_f0_quadrature(f0_bounds, n_nodes=256):
    """Return f0 nodes and normalized weights for a log-uniform f0 prior."""
    bounds = np.asarray(f0_bounds, dtype=float)
    if bounds.shape != (2,) or not np.isfinite(bounds).all():
        raise ValueError("f0_bounds must contain two finite values")
    low, high = map(float, bounds)
    ref_low, ref_high = F0_REFERENCE_BOUNDS
    if not 0.0 < ref_low <= low < high <= ref_high:
        raise ValueError("f0_bounds must lie inside the reference f0 envelope")
    if isinstance(n_nodes, (bool, np.bool_)) or int(n_nodes) != n_nodes or n_nodes < 2:
        raise ValueError("n_nodes must be an integer of at least 2")

    log_span = np.log(ref_high) - np.log(ref_low)
    u_low = (np.log(low) - np.log(ref_low)) / log_span
    u_high = (np.log(high) - np.log(ref_low)) / log_span
    nodes, weights = np.polynomial.legendre.leggauss(int(n_nodes))
    u = u_low + (nodes + 1.0) * (u_high - u_low) / 2.0
    f0 = np.exp(np.log(ref_low) + u * log_span)
    return f0, weights / 2.0


def marginal_signal_log_density(
    observations: SignalObservations,
    parameters: SignalModelParameters,
    *,
    f0_bounds,
    n_nodes=256,
):
    """Integrate p(molecular report, MRI score | f0) over the declared f0 prior.

    The returned vector has one log density per paired observation. The f0 prior is
    log-uniform over ``f0_bounds`` and is shared by all rows. The model parameters are
    scenario assumptions, not fitted values.
    """
    if not all(
        hasattr(observations, name)
        for name in ("molecular_report", "mri_score")
    ):
        raise TypeError("observations must provide molecular_report and mri_score")
    reports = np.asarray(observations.molecular_report)
    scores = np.asarray(observations.mri_score, dtype=float)
    if reports.ndim != 1 or scores.shape != reports.shape or reports.size < 1:
        raise ValueError("signal arrays must be nonempty one-dimensional arrays of equal size")
    f0_nodes, weights = log_uniform_f0_quadrature(f0_bounds, n_nodes=n_nodes)
    particles = np.broadcast_to(f0_nodes, (reports.size, f0_nodes.size))
    conditional = integrated_signal_log_likelihood(
        particles,
        SignalObservations(reports, scores),
        parameters,
    )
    return logsumexp(conditional + np.log(weights)[None, :], axis=1)


def compare_marginal_signal_models(
    parameters_p: SignalModelParameters,
    parameters_q: SignalModelParameters,
    *,
    f0_bounds,
    n_nodes=256,
    n_mri=1601,
):
    """Compare two marginal signal models by normalization, KL, and max density gap.

    The integration uses a fixed MRI-score grid from -8 to 8 standard units. This is a
    deterministic model-identifiability diagnostic, not a finite-sample power analysis.
    """
    if isinstance(n_mri, (bool, np.bool_)) or int(n_mri) != n_mri or n_mri < 2:
        raise ValueError("n_mri must be an integer of at least 2")
    scores = np.linspace(-8.0, 8.0, int(n_mri))
    p_densities, q_densities = [], []
    for report in (False, True):
        observations = SignalObservations(
            molecular_report=np.full(scores.size, report, dtype=bool),
            mri_score=scores,
        )
        p_densities.append(
            np.exp(
                marginal_signal_log_density(
                    observations,
                    parameters_p,
                    f0_bounds=f0_bounds,
                    n_nodes=n_nodes,
                )
            )
        )
        q_densities.append(
            np.exp(
                marginal_signal_log_density(
                    observations,
                    parameters_q,
                    f0_bounds=f0_bounds,
                    n_nodes=n_nodes,
                )
            )
        )
    p_densities = np.asarray(p_densities)
    q_densities = np.asarray(q_densities)
    epsilon = np.finfo(float).tiny
    kl_pq = float(
        sum(
            np.trapezoid(
                p * (np.log(np.maximum(p, epsilon)) - np.log(np.maximum(q, epsilon))),
                scores,
            )
            for p, q in zip(p_densities, q_densities)
        )
    )
    kl_qp = float(
        sum(
            np.trapezoid(
                q * (np.log(np.maximum(q, epsilon)) - np.log(np.maximum(p, epsilon))),
                scores,
            )
            for p, q in zip(p_densities, q_densities)
        )
    )
    return {
        "p_mass": float(np.trapezoid(p_densities.sum(axis=0), scores)),
        "q_mass": float(np.trapezoid(q_densities.sum(axis=0), scores)),
        "kl_pq": kl_pq,
        "kl_qp": kl_qp,
        "max_abs_density_difference": float(
            np.max(np.abs(p_densities - q_densities))
        ),
        "n_f0_nodes": int(n_nodes),
        "n_mri_grid": int(n_mri),
    }


def main():
    """Print a deterministic comparison for one equal-signature parameter pair."""
    import json

    first = SignalModelParameters(q=0.8, a=0.5, b=0.4, d=0.6)
    q, rho = standard_normal_marginal_signature(first)
    a_alt, b_alt = 0.8, 0.5
    d_alt = (rho - a_alt * b_alt) / (
        np.sqrt(1.0 - a_alt * a_alt) * np.sqrt(1.0 - b_alt * b_alt)
    )
    second = SignalModelParameters(q=q, a=a_alt, b=b_alt, d=d_alt)
    f0_ranges = {
        "calibration_low_shift": (
            0.001 * np.exp(-0.7), 0.1 * np.exp(-0.7)
        ),
        "default": (0.001, 0.1),
        "calibration_high_shift": (0.001 * np.exp(2.3), 0.9),
    }
    result = {
        "scope": "deterministic synthetic likelihood diagnostic; not policy evaluation",
        "assumptions": {
            "q": q,
            "f0_prior": "log-uniform on each listed range",
            "f0_quadrature_nodes": 256,
            "mri_grid": [-8.0, 8.0],
            "mri_grid_points": 1601,
        },
        "normal_score_signature": {"q": q, "rho": rho},
        "parameter_pair": {
            "first": {"a": first.a, "b": first.b, "d": first.d},
            "second": {"a": second.a, "b": second.b, "d": second.d},
        },
        "by_f0_range": {
            name: {
                "f0_bounds": [float(low), float(high)],
                **compare_marginal_signal_models(
                    first,
                    second,
                    f0_bounds=(low, high),
                    n_nodes=256,
                    n_mri=1601,
                ),
            }
            for name, (low, high) in f0_ranges.items()
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True))


__all__ = [
    "compare_marginal_signal_models",
    "log_uniform_f0_quadrature",
    "marginal_signal_log_density",
]


if __name__ == "__main__":
    main()
