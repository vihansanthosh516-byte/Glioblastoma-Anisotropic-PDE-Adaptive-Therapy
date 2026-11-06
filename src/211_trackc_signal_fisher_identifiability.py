"""Local Fisher-information diagnostic for assumed Track C signal parameters.

This is deterministic and model-based. It does not estimate parameters from patients
or evaluate treatment policies. See the frozen protocol in the novelty update.
"""
from __future__ import annotations

import json
import importlib.util
import sys
from pathlib import Path

import numpy as np

from trackc_integrated_signal_model import (
    F0_REFERENCE_BOUNDS,
    SignalModelParameters,
    SignalObservations,
)
_DEPENDENCE_PATH = Path(__file__).with_name("210_trackc_dependence_identifiability.py")
_DEPENDENCE_SPEC = importlib.util.spec_from_file_location(
    "trackc_dependence_identifiability", _DEPENDENCE_PATH
)
if _DEPENDENCE_SPEC is None or _DEPENDENCE_SPEC.loader is None:
    raise ImportError(f"cannot load dependence likelihood module: {_DEPENDENCE_PATH}")
_DEPENDENCE_MODULE = importlib.util.module_from_spec(_DEPENDENCE_SPEC)
sys.modules[_DEPENDENCE_SPEC.name] = _DEPENDENCE_MODULE
_DEPENDENCE_SPEC.loader.exec_module(_DEPENDENCE_MODULE)
marginal_signal_log_density = _DEPENDENCE_MODULE.marginal_signal_log_density


PARAMETER_NAMES = ("q", "a", "b", "d")
NOMINAL_PARAMETERS = (0.8, 0.5, 0.4, 0.6)
CALIBRATION_F0_RANGES = {
    "calibration_low_shift": (0.001 * float(np.exp(-0.7)), 0.1 * float(np.exp(-0.7))),
    "default": (0.001, 0.1),
    "calibration_high_shift": (0.001 * float(np.exp(2.3)), 0.9),
}


def _signal_density(parameters, f0_bounds, *, n_f0_nodes, mri_grid):
    densities = []
    for report in (False, True):
        observations = SignalObservations(
            molecular_report=np.full(mri_grid.size, report, dtype=bool),
            mri_score=mri_grid,
        )
        densities.append(
            np.exp(
                marginal_signal_log_density(
                    observations,
                    parameters,
                    f0_bounds=f0_bounds,
                    n_nodes=n_f0_nodes,
                )
            )
        )
    return np.asarray(densities)


def expected_fisher_information(
    f0_bounds,
    *,
    parameters=NOMINAL_PARAMETERS,
    n_f0_nodes=256,
    n_mri=1601,
    finite_difference_step=1e-4,
    rank_relative_tolerance=1e-6,
):
    """Compute per-pair Fisher information for (q, a, b, d).

    The likelihood marginalizes over a log-uniform f0 prior. MRI-score density is
    integrated over [-8, 8] with a trapezoid rule. The numerical rank threshold is
    fixed by the predeclared protocol.
    """
    theta = np.asarray(parameters, dtype=float)
    if theta.shape != (4,) or not np.isfinite(theta).all():
        raise ValueError("parameters must be four finite values (q, a, b, d)")
    if not 0.0 < theta[0] < 1.0 or np.any(np.abs(theta[1:]) >= 1.0):
        raise ValueError("parameters must be interior to q in (0,1), a/b/d in (-1,1)")
    if not np.isfinite(finite_difference_step) or finite_difference_step <= 0.0:
        raise ValueError("finite_difference_step must be finite and positive")
    if not 0.0 < rank_relative_tolerance < 1.0:
        raise ValueError("rank_relative_tolerance must lie in (0,1)")
    if isinstance(n_mri, (bool, np.bool_)) or int(n_mri) != n_mri or n_mri < 3:
        raise ValueError("n_mri must be an integer of at least 3")

    mri_grid = np.linspace(-8.0, 8.0, int(n_mri))
    base_parameters = SignalModelParameters(*map(float, theta))
    density = _signal_density(
        base_parameters, f0_bounds, n_f0_nodes=n_f0_nodes, mri_grid=mri_grid
    )
    if not np.isfinite(density).all() or np.any(density <= 0.0):
        raise ValueError("marginal signal density must be finite and positive on grid")

    derivatives = np.empty((4, 2, mri_grid.size), dtype=float)
    for index in range(4):
        plus = theta.copy()
        minus = theta.copy()
        plus[index] += finite_difference_step
        minus[index] -= finite_difference_step
        if index == 0 and not (0.0 < minus[index] < plus[index] < 1.0):
            raise ValueError("finite-difference step crosses q boundary")
        if index > 0 and np.any(np.abs([minus[index], plus[index]]) >= 1.0):
            raise ValueError("finite-difference step crosses coefficient boundary")
        p_plus = _signal_density(
            SignalModelParameters(*map(float, plus)),
            f0_bounds,
            n_f0_nodes=n_f0_nodes,
            mri_grid=mri_grid,
        )
        p_minus = _signal_density(
            SignalModelParameters(*map(float, minus)),
            f0_bounds,
            n_f0_nodes=n_f0_nodes,
            mri_grid=mri_grid,
        )
        derivatives[index] = (p_plus - p_minus) / (2.0 * finite_difference_step)

    fisher = np.empty((4, 4), dtype=float)
    for i in range(4):
        for j in range(4):
            integrand = np.sum(derivatives[i] * derivatives[j] / density, axis=0)
            fisher[i, j] = np.trapezoid(integrand, mri_grid)
    fisher = (fisher + fisher.T) / 2.0
    eigenvalues = np.linalg.eigvalsh(fisher)
    largest = float(max(eigenvalues[-1], 0.0))
    effective_rank = int(np.count_nonzero(eigenvalues > rank_relative_tolerance * largest))
    condition_number = (
        float(eigenvalues[-1] / eigenvalues[0])
        if effective_rank == 4 and eigenvalues[0] > 0.0
        else None
    )
    d_standard_error_bound = (
        float(np.sqrt(np.linalg.inv(fisher)[3, 3]))
        if effective_rank == 4
        else None
    )
    mass = float(np.trapezoid(density.sum(axis=0), mri_grid))
    return {
        "f0_bounds": [float(f0_bounds[0]), float(f0_bounds[1])],
        "parameter_names": list(PARAMETER_NAMES),
        "parameters": [float(x) for x in theta],
        "f0_nodes": int(n_f0_nodes),
        "mri_grid": [-8.0, 8.0],
        "mri_grid_points": int(n_mri),
        "finite_difference_step": float(finite_difference_step),
        "rank_relative_tolerance": float(rank_relative_tolerance),
        "density_mass": mass,
        "fisher_information_per_pair": fisher.tolist(),
        "eigenvalues": [float(x) for x in eigenvalues],
        "effective_rank": effective_rank,
        "condition_number": condition_number,
        "d_local_cramer_rao_sd_bound_per_sqrt_n": d_standard_error_bound,
    }


def build_diagnostic():
    """Evaluate reference-normal control and three fixed calibration ranges."""
    ranges = {"reference_normal_control": F0_REFERENCE_BOUNDS, **CALIBRATION_F0_RANGES}
    return {
        "scope": "deterministic local Fisher-information diagnostic; no patient data or policy outcomes",
        "limitations": [
            "one assumed parameter point",
            "local information is not finite-sample recovery or power",
            "molecular-to-resistance mapping is assumed in simulation",
        ],
        "by_f0_range": {
            name: expected_fisher_information(bounds)
            for name, bounds in ranges.items()
        },
    }


def main():
    print(json.dumps(build_diagnostic(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["build_diagnostic", "expected_fisher_information"]
