"""Predeclared signal-strength grid stress test for Track C identifiability.

This deterministic diagnostic evaluates every frozen grid point. It does not select
parameters, estimate patient biology, or compare treatment outcomes.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np

from trackc_integrated_signal_model import F0_REFERENCE_BOUNDS

_FISHER_PATH = Path(__file__).with_name("211_trackc_signal_fisher_identifiability.py")
_FISHER_SPEC = importlib.util.spec_from_file_location(
    "trackc_signal_fisher_identifiability", _FISHER_PATH
)
if _FISHER_SPEC is None or _FISHER_SPEC.loader is None:
    raise ImportError(f"cannot load Fisher-information module: {_FISHER_PATH}")
_FISHER_MODULE = importlib.util.module_from_spec(_FISHER_SPEC)
sys.modules[_FISHER_SPEC.name] = _FISHER_MODULE
_FISHER_SPEC.loader.exec_module(_FISHER_MODULE)
expected_fisher_information = _FISHER_MODULE.expected_fisher_information


Q_GRID = (0.5, 0.7, 0.9)
COEFFICIENT_GRID = (-0.5, 0.0, 0.5)
F0_RANGES = {
    "reference_normal_control": F0_REFERENCE_BOUNDS,
    **_FISHER_MODULE.CALIBRATION_F0_RANGES,
}
N_F0_NODES = 256
N_MRI_POINTS = 1601
FINITE_DIFFERENCE_STEP = 1e-4
RANK_RELATIVE_TOLERANCE = 1e-6


def parameter_grid():
    """Return all 81 frozen (q, a, b, d) combinations in stable order."""
    return [
        (float(q), float(a), float(b), float(d))
        for q, a, b, d in product(Q_GRID, COEFFICIENT_GRID, COEFFICIENT_GRID, COEFFICIENT_GRID)
    ]


def build_diagnostic(*, n_f0_nodes=N_F0_NODES, n_mri=N_MRI_POINTS):
    """Evaluate every frozen parameter point under all four f0 distributions."""
    results = []
    for range_name, bounds in F0_RANGES.items():
        for parameters in parameter_grid():
            summary = expected_fisher_information(
                bounds,
                parameters=parameters,
                n_f0_nodes=n_f0_nodes,
                n_mri=n_mri,
                finite_difference_step=FINITE_DIFFERENCE_STEP,
                rank_relative_tolerance=RANK_RELATIVE_TOLERANCE,
            )
            if abs(summary["density_mass"] - 1.0) > 1e-5:
                raise ArithmeticError(
                    f"marginal density does not normalize for {range_name}, {parameters}"
                )
            results.append(
                {
                    "f0_range": range_name,
                    "parameters": list(parameters),
                    "effective_rank": summary["effective_rank"],
                    "condition_number": summary["condition_number"],
                    "d_local_cramer_rao_sd_bound_per_sqrt_n": summary[
                        "d_local_cramer_rao_sd_bound_per_sqrt_n"
                    ],
                    "eigenvalues": summary["eigenvalues"],
                    "density_mass": summary["density_mass"],
                }
            )

    by_range = {}
    for range_name in F0_RANGES:
        rows = [row for row in results if row["f0_range"] == range_name]
        rank_counts = {
            str(rank): sum(row["effective_rank"] == rank for row in rows)
            for rank in range(5)
        }
        full_rank_bounds = [
            row["d_local_cramer_rao_sd_bound_per_sqrt_n"]
            for row in rows
            if row["d_local_cramer_rao_sd_bound_per_sqrt_n"] is not None
        ]
        by_range[range_name] = {
            "n_points": len(rows),
            "effective_rank_counts": rank_counts,
            "full_rank_d_bound_min_median_max": (
                [
                    float(np.min(full_rank_bounds)),
                    float(np.median(full_rank_bounds)),
                    float(np.max(full_rank_bounds)),
                ]
                if full_rank_bounds
                else None
            ),
            "q_0_5_all_d_unidentified": all(
                row["d_local_cramer_rao_sd_bound_per_sqrt_n"] is None
                for row in rows
                if row["parameters"][0] == 0.5
            ),
        }

    return {
        "scope": "deterministic synthetic local-information grid; no patient data or treatment outcomes",
        "protocol": {
            "q_grid": list(Q_GRID),
            "a_b_d_grid": list(COEFFICIENT_GRID),
            "parameter_points": len(parameter_grid()),
            "f0_ranges": {
                name: [float(bounds[0]), float(bounds[1])]
                for name, bounds in F0_RANGES.items()
            },
            "f0_quadrature_nodes": int(n_f0_nodes),
            "mri_grid": [-8.0, 8.0],
            "mri_grid_points": int(n_mri),
            "finite_difference_step": FINITE_DIFFERENCE_STEP,
            "rank_relative_tolerance": RANK_RELATIVE_TOLERANCE,
        },
        "summary_by_f0_range": by_range,
        "all_points": results,
    }


def main():
    print(json.dumps(build_diagnostic(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = ["build_diagnostic", "parameter_grid"]
