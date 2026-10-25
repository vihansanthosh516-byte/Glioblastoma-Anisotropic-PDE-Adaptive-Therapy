"""Synthetic identifiability check for Track B's cross-fitted boundary variogram estimator.

This isolates covariance estimation from PDE fitting and thresholded-mask generation. Each synthetic patient
contributes one independent Gaussian residual field drawn exactly from the exponential-plus-nugget model that
the estimator assumes. It uses no patient data and reports only aggregate recovery and boundary-selection rates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.fk_boundary_calibration import fit_cross_fitted_boundary_covariance
from src.fk_spatial_uncertainty import deterministic_surface_anchors, exponential_covariance_matrix


DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "track_b_calibration_identifiability_20261024.json"
TRUTH = {"sigma_b_mm": 0.6, "tau_b_mm": 0.2, "ell_mm": 4.0}
CURRENT_CANDIDATES = (5.0, 10.0, 20.0)
WIDE_CANDIDATES = (1.0, 2.0, 4.0, 5.0, 10.0, 20.0, 40.0, 80.0)


def _surface_anchors(n_anchors: int) -> np.ndarray:
    shape = (32, 32, 32)
    coords = np.ogrid[tuple(slice(0, size) for size in shape)]
    center = tuple((size - 1) / 2.0 for size in shape)
    sphere = sum((axis - center_i) ** 2 for axis, center_i in zip(coords, center)) <= 8.0**2
    return deterministic_surface_anchors(sphere, spacing_mm=1.0, n_anchors=n_anchors)


def _fit_summary(estimates: list[dict], candidates: tuple[float, ...]) -> dict:
    lengths = np.asarray([row["ell_mm"] for row in estimates], dtype=float)
    sigmas = np.asarray([row["sigma_b_mm"] for row in estimates], dtype=float)
    nuggets = np.asarray([row["tau_b_mm"] for row in estimates], dtype=float)
    return {
        "candidate_ell_mm": list(candidates),
        "true_ell_in_candidate_grid": bool(TRUTH["ell_mm"] in candidates),
        "selected_at_lower_bound": int(np.count_nonzero(lengths == min(candidates))),
        "selected_at_upper_bound": int(np.count_nonzero(lengths == max(candidates))),
        "selected_ell_median_mm": float(np.median(lengths)),
        "selected_ell_q10_q90_mm": [float(v) for v in np.quantile(lengths, [0.1, 0.9])],
        "selected_ell_median_absolute_error_mm": float(np.median(np.abs(lengths - TRUTH["ell_mm"]))),
        "sigma_b_median_mm": float(np.median(sigmas)),
        "sigma_b_median_absolute_error_mm": float(np.median(np.abs(sigmas - TRUTH["sigma_b_mm"]))),
        "tau_b_median_mm": float(np.median(nuggets)),
        "tau_b_median_absolute_error_mm": float(np.median(np.abs(nuggets - TRUTH["tau_b_mm"]))),
        "exact_true_ell_selection_count": int(np.count_nonzero(lengths == TRUTH["ell_mm"])),
    }


def run_identifiability_check(*, replicates: int = 30, patient_counts=(8, 32),
                              anchor_counts=(24, 256), seed: int = 20261029) -> dict:
    if isinstance(replicates, bool) or int(replicates) != replicates or replicates < 2:
        raise ValueError("replicates must be an integer of at least two")
    if isinstance(seed, bool) or int(seed) != seed:
        raise ValueError("seed must be an integer")
    patient_counts = tuple(int(value) for value in patient_counts)
    anchor_counts = tuple(int(value) for value in anchor_counts)
    if (not patient_counts or any(value < 3 for value in patient_counts)
            or not anchor_counts or any(value < 3 for value in anchor_counts)):
        raise ValueError("each scenario needs at least three patients and three anchors")

    rng = np.random.default_rng(int(seed))
    scenarios = []
    for n_patients in patient_counts:
        for n_anchors in anchor_counts:
            anchors = _surface_anchors(n_anchors)
            covariance = exponential_covariance_matrix(anchors, **TRUTH)
            factor = np.linalg.cholesky(covariance)
            current_estimates = []
            wide_estimates = []
            for replicate in range(int(replicates)):
                residual_matrix = factor @ rng.standard_normal((len(anchors), n_patients))
                cases = []
                for patient_index in range(n_patients):
                    patient_id = f"synthetic-{patient_index:03d}"
                    training_ids = [f"synthetic-{other:03d}" for other in range(n_patients)
                                    if other != patient_index]
                    cases.append({
                        "patient_id": patient_id,
                        "origin_id": f"replicate-{replicate:03d}",
                        "training_patient_ids": training_ids,
                        "residuals_mm": residual_matrix[:, patient_index],
                        "anchor_points_mm": anchors,
                    })
                current_estimates.append(fit_cross_fitted_boundary_covariance(
                    cases, ell_candidates_mm=CURRENT_CANDIDATES
                ))
                wide_estimates.append(fit_cross_fitted_boundary_covariance(
                    cases, ell_candidates_mm=WIDE_CANDIDATES
                ))
            scenarios.append({
                "n_calibration_patients": n_patients,
                "n_surface_anchors_per_patient": int(len(anchors)),
                "replicates": int(replicates),
                "current_grid": _fit_summary(current_estimates, CURRENT_CANDIDATES),
                "wide_grid": _fit_summary(wide_estimates, WIDE_CANDIDATES),
            })

    return {
        "status": "synthetic_estimator_diagnostic",
        "seed": int(seed),
        "truth": dict(TRUTH),
        "covariance_model": "sigma_b^2 exp(-distance/ell) + tau_b^2 I",
        "anchor_geometry": "deterministically spaced voxels on a voxelized 8 mm radius sphere at 1 mm spacing",
        "scenarios": scenarios,
        "interpretation_limit": (
            "This isolates the estimator under its exactly correct Gaussian covariance model. It does not include "
            "PDE fitting, model discrepancy, binary-mask thresholding, patient scans, or prove that the model fits "
            "real segmentation residuals."
        ),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    report = run_identifiability_check(replicates=args.replicates, seed=args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
