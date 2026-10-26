"""Synthetic split-conformal reliability check for Track B boundary forecasts.

The generator and predictor use no patient data. The calibration and test partitions are patient-disjoint.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy import ndimage

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src import fk_full_hybrid_gate as gate
from src.fk_hybrid_forecast import hybrid_posterior_predictive_forecast
from src.fk_patient_cv import ForecastFailure, PatientHistory, PatientSeries
from src.fk_patient_conformal import (
    split_conformal_radius,
    symmetric_boundary_hausdorff_mm,
    target_boundary_coverage,
)


FORECAST_CONFIG = {"n_coarse": 3, "n_refine": 3, "n_draws": 2, "n_anchors": 24}
ARMS = {
    "oracle_covariance": {
        "covariance": {"sigma_b_mm": 0.6, "tau_b_mm": 0.2, "ell_mm": 4.0},
        "shift_voxels_xyz": (0, 0, 0),
    },
    "endpoint_covariance": {
        "covariance": {"sigma_b_mm": 3.4, "tau_b_mm": 0.3, "ell_mm": 80.0},
        "shift_voxels_xyz": (0, 0, 0),
    },
    "endpoint_covariance_plus_2mm_bias": {
        "covariance": {"sigma_b_mm": 3.4, "tau_b_mm": 0.3, "ell_mm": 80.0},
        "shift_voxels_xyz": (0, 0, 1),
    },
    "endpoint_covariance_plus_2mm_bias_train_corrected": {
        "covariance": {"sigma_b_mm": 3.4, "tau_b_mm": 0.3, "ell_mm": 80.0},
        "shift_voxels_xyz": (0, 0, 1),
        "training_centroid_correction": True,
    },
}


def _iid_mixture_manifest(*, seed: int, n_patients: int = 60,
                          factor_design: str = "independent") -> list[dict]:
    """Sample conditions, scan-gap order, and anisotropy axis from an explicit design."""
    if factor_design not in {"independent", "legacy_coupled"}:
        raise ValueError("factor_design must be 'independent' or 'legacy_coupled'")
    rng = np.random.default_rng(seed)
    rows = []
    for index in range(n_patients):
        diffusion = str(rng.choice(gate.DIFFUSION_TRUTHS))
        volume_trend = str(rng.choice(gate.VOLUME_TRENDS))
        replicate = int(rng.integers(0, 15))
        if factor_design == "independent":
            factor_rng = np.random.default_rng(
                np.random.SeedSequence([int(seed), index, 791])
            )
            gap_pattern_index = int(factor_rng.integers(0, 3))
            anisotropy_axis_index = int(factor_rng.integers(0, 3))
        else:
            gap_pattern_index = replicate % 3
            anisotropy_axis_index = replicate % 3
        case_seed = int(np.random.SeedSequence([int(seed), index, 91]).generate_state(1)[0])
        rows.append({
            "case_id": f"mix-{seed}-{index:03d}",
            "diffusion_truth": diffusion,
            "volume_trend": volume_trend,
            "replicate": replicate,
            "gap_pattern_index": gap_pattern_index,
            "anisotropy_axis_index": anisotropy_axis_index,
            "seed": case_seed,
        })
    return rows


def _random_split(manifest: list[dict], *, seed: int) -> dict[str, list[str]]:
    """Randomly split iid synthetic patients so calibration and test scores remain exchangeable."""
    if len(manifest) != 60:
        raise ValueError("this experiment requires exactly 60 synthetic patients")
    order = np.random.default_rng(seed).permutation(len(manifest))
    ids = [manifest[index]["case_id"] for index in order]
    return {"train": ids[:20], "calibration": ids[20:40], "test": ids[40:]}


def _first_third_visit(record: PatientSeries) -> PatientSeries:
    return PatientSeries(record.patient_id, np.asarray(record.days[:3], dtype=float), tuple(record.masks[:3]))


def _predict_boundary(record: PatientSeries, nuisance: dict, *, shift_voxels_xyz) -> tuple[np.ndarray, bool]:
    history = PatientHistory(
        patient_id=record.patient_id,
        target_index=2,
        forecast_day=float(record.days[2]),
        days=np.asarray(record.days[:2], dtype=float),
        masks=tuple(record.masks[:2]),
    )
    try:
        prediction = gate._forecast_callback(
            history, nuisance,
            forecast_function=hybrid_posterior_predictive_forecast,
            forecast_config=FORECAST_CONFIG,
        )
        draws = np.asarray(prediction["forecast"]["mask_draws"], dtype=bool)
        point_mask = draws.mean(axis=0) >= 0.5
        fallback = False
    except ForecastFailure:
        point_mask = np.asarray(history.masks[-1], dtype=bool).copy()
        fallback = True
    if shift_voxels_xyz != (0, 0, 0):
        point_mask = ndimage.shift(
            point_mask.astype(np.uint8), shift=shift_voxels_xyz,
            order=0, mode="constant", cval=0,
        ).astype(bool)
    if not point_mask.any() or point_mask.all():
        point_mask = np.asarray(history.masks[-1], dtype=bool).copy()
        fallback = True
    return point_mask, fallback


def _mask_centroid_mm(mask) -> np.ndarray:
    coords = np.argwhere(np.asarray(mask, dtype=bool))
    if not len(coords):
        raise ValueError("cannot find centroid of an empty mask")
    return coords.mean(axis=0) * gate.SPACING_MM


def _estimate_training_centroid_correction(training_records, *, covariance, imposed_shift, seed: int) -> dict:
    """Estimate a global translation from two-fold held-out training patients only."""
    records = tuple(training_records)
    if len(records) < 4:
        raise ValueError("training-only centroid correction requires at least four patients")
    order = np.random.default_rng(seed).permutation(len(records))
    validation_folds = np.array_split(order, 2)
    offsets_voxels = []
    fallback_count = 0
    for validation_indices in validation_folds:
        validation_set = set(int(index) for index in validation_indices)
        fit_records = [record for index, record in enumerate(records) if index not in validation_set]
        validation_records = [records[index] for index in validation_indices]
        nuisance = gate.fit_training_nuisance(fit_records)
        nuisance["boundary_covariance"] = dict(covariance)
        for record in validation_records:
            predicted, fallback = _predict_boundary(
                record, nuisance, shift_voxels_xyz=imposed_shift
            )
            fallback_count += int(fallback)
            target = np.asarray(record.masks[2], dtype=bool)
            offsets_voxels.append(
                (_mask_centroid_mm(target) - _mask_centroid_mm(predicted)) / gate.SPACING_MM
            )
    correction = np.rint(np.median(np.stack(offsets_voxels), axis=0)).astype(int)
    return {
        "correction_voxels_xyz": tuple(int(value) for value in correction),
        "n_cross_fitted_training_patients": len(offsets_voxels),
        "n_forecast_fallbacks": int(fallback_count),
    }


def _run_cohort(*, seed: int, alphas=(0.1, 0.05), factor_design="independent") -> dict:
    manifest = _iid_mixture_manifest(seed=seed, factor_design=factor_design)
    split = _random_split(manifest, seed=seed + 3)
    generated = {record.patient_id: (_first_third_visit(record), truth)
                 for row in manifest for record, truth in [gate.generate_synthetic_series(row)]}
    training = [generated[patient_id][0] for patient_id in split["train"]]
    nuisance_base = gate.fit_training_nuisance(training)
    rows_by_id = {row["case_id"]: row for row in manifest}
    condition_counts = {}
    for split_name, patient_ids in split.items():
        counts: dict[str, int] = {}
        for patient_id in patient_ids:
            row = rows_by_id[patient_id]
            condition = f"{row['diffusion_truth']}|{row['volume_trend']}"
            counts[condition] = counts.get(condition, 0) + 1
        condition_counts[split_name] = counts
    factor_counts = {}
    for split_name, patient_ids in split.items():
        rows = [rows_by_id[patient_id] for patient_id in patient_ids]
        factor_counts[split_name] = {
            f"gap{gap}_axis{axis}": sum(
                row["gap_pattern_index"] == gap and row["anisotropy_axis_index"] == axis
                for row in rows
            )
            for gap in range(3) for axis in range(3)
        }
    result = {"seed": int(seed), "factor_design": factor_design,
              "partition_sizes": {key: len(value) for key, value in split.items()},
              "condition_counts_by_split": condition_counts,
              "gap_axis_counts_by_split": factor_counts, "arms": {}}

    for arm_name, arm in ARMS.items():
        nuisance = dict(nuisance_base)
        nuisance["boundary_covariance"] = dict(arm["covariance"])
        correction = {"correction_voxels_xyz": (0, 0, 0),
                      "n_cross_fitted_training_patients": 0,
                      "n_forecast_fallbacks": 0}
        if arm.get("training_centroid_correction", False):
            correction = _estimate_training_centroid_correction(
                training, covariance=arm["covariance"],
                imposed_shift=arm["shift_voxels_xyz"], seed=seed + 17,
            )
        effective_shift = tuple(
            int(base + adjustment)
            for base, adjustment in zip(arm["shift_voxels_xyz"], correction["correction_voxels_xyz"])
        )
        cal_rows = []
        test_rows = []
        fallback_count = int(correction["n_forecast_fallbacks"])
        for patient_id in split["calibration"]:
            record, _ = generated[patient_id]
            predicted, fallback = _predict_boundary(
                record, nuisance, shift_voxels_xyz=effective_shift
            )
            fallback_count += int(fallback)
            target = np.asarray(record.masks[2], dtype=bool)
            cal_rows.append({"score_mm": symmetric_boundary_hausdorff_mm(
                predicted, target, spacing_mm=gate.SPACING_MM
            )})
        for patient_id in split["test"]:
            record, _ = generated[patient_id]
            predicted, fallback = _predict_boundary(
                record, nuisance, shift_voxels_xyz=effective_shift
            )
            fallback_count += int(fallback)
            target = np.asarray(record.masks[2], dtype=bool)
            test_rows.append({
                "score_mm": symmetric_boundary_hausdorff_mm(
                    predicted, target, spacing_mm=gate.SPACING_MM
                ),
                "predicted": predicted,
                "target": target,
            })

        arm_result = {"n_calibration_patients": len(cal_rows), "n_test_patients": len(test_rows),
                      "training_centroid_correction_voxels_xyz": list(correction["correction_voxels_xyz"]),
                      "n_cross_fitted_training_patients_for_correction": int(
                          correction["n_cross_fitted_training_patients"]
                      ),
                      "n_forecast_fallbacks": fallback_count,
                      "median_raw_hausdorff_mm": float(np.median([row["score_mm"] for row in test_rows])),
                      "coverage": {}}
        scores = [row["score_mm"] for row in cal_rows]
        for alpha in alphas:
            quantile = split_conformal_radius(scores, alpha=alpha)
            coverage_rows = [target_boundary_coverage(
                row["predicted"], row["target"], radius_mm=quantile["radius_mm"],
                spacing_mm=gate.SPACING_MM,
            ) for row in test_rows]
            level = int(round(100 * (1.0 - alpha)))
            arm_result["coverage"][str(level)] = {
                "calibration_rank": quantile["rank"],
                "radius_mm": float(quantile["radius_mm"]),
                "finite_radius": bool(quantile["finite"]),
                "whole_boundary_patient_count": int(sum(
                    row["whole_boundary_covered"] for row in coverage_rows
                )),
                "whole_boundary_coverage": float(np.mean([
                    row["whole_boundary_covered"] for row in coverage_rows
                ])),
                "mean_target_surface_fraction_covered": float(np.mean([
                    row["target_surface_fraction_covered"] for row in coverage_rows
                ])),
                "median_band_voxel_fraction": float(np.median([
                    row["band_voxel_fraction"] for row in coverage_rows
                ])),
            }
        result["arms"][arm_name] = arm_result
    return result


def run_conformal_synthetic_gate(*, replicates: int = 5, seed: int = 20261030,
                                 factor_design: str = "independent") -> dict:
    if isinstance(replicates, bool) or int(replicates) != replicates or replicates < 1:
        raise ValueError("replicates must be a positive integer")
    if factor_design not in {"independent", "legacy_coupled"}:
        raise ValueError("factor_design must be 'independent' or 'legacy_coupled'")
    cohorts = []
    for index in range(int(replicates)):
        cohorts.append(_run_cohort(seed=int(seed) + index, factor_design=factor_design))
        print(f"completed cohort {index + 1}/{int(replicates)}", flush=True)
    aggregate = {}
    for arm_name in ARMS:
        aggregate[arm_name] = {}
        for level in (90, 95):
            rows = [cohort["arms"][arm_name]["coverage"][str(level)] for cohort in cohorts]
            scores = [cohort["arms"][arm_name] for cohort in cohorts]
            rates = [row["whole_boundary_coverage"] for row in rows]
            radii = [row["radius_mm"] for row in rows]
            aggregate[arm_name][str(level)] = {
                "n_cohorts": len(cohorts),
                "mean_cohort_whole_boundary_coverage": float(np.mean(rates)),
                "minimum_cohort_whole_boundary_coverage": float(np.min(rates)),
                "maximum_cohort_whole_boundary_coverage": float(np.max(rates)),
                "median_calibrated_radius_mm": float(np.median(radii)),
                "median_band_voxel_fraction": float(np.median([
                    row["median_band_voxel_fraction"] for row in rows
                ])),
                "total_whole_boundary_covered_test_patients": int(sum(
                    row["whole_boundary_patient_count"] for row in rows
                )),
                "total_test_patients": int(sum(cohort["arms"][arm_name]["n_test_patients"]
                                                for cohort in cohorts)),
                "total_forecast_fallbacks": int(sum(score["n_forecast_fallbacks"] for score in scores)),
            }
    return {
        "status": "synthetic_split_conformal_development_pilot",
        "seed": int(seed),
        "cohort_replicates": int(replicates),
        "factor_design": factor_design,
        "patient_conditions": (
            "independent uniform draws from four synthetic conditions; scan-gap order and anisotropic axis are "
            + ("sampled independently per patient; patients split at random" if factor_design == "independent"
               else "both derived from replicate modulo three as in the historical synthetic runs")
        ),
        "split_per_cohort": {"training_patients": 20, "calibration_patients": 20, "test_patients": 20},
        "forecast_origin": "first two masks predict third mask; one target per patient",
        "forecast_config": dict(FORECAST_CONFIG),
        "runtime_note": "reduced PDE search grid and two draws for a CPU development pilot; not a final accuracy benchmark",
        "covariance_arms": {name: {"covariance": arm["covariance"],
                                    "systematic_shift_voxels_xyz": list(arm["shift_voxels_xyz"]),
                                    "training_centroid_correction": bool(
                                        arm.get("training_centroid_correction", False)
                                    )}
                            for name, arm in ARMS.items()},
        "score": "symmetric Hausdorff distance between predicted and target six-connected boundary surfaces, in mm",
        "conformal_method": "patient-level split conformal; quantile rank ceil((n_calibration+1)*(1-alpha))",
        "aggregate": aggregate,
        "cohorts": cohorts,
        "interpretation_limit": (
            "synthetic development only; the calibration guarantee is marginal for exchangeable patients when the "
            "predictor is frozen before calibration. The generator is not an MRI noise model, and the reported "
            "band width has no established clinical acceptability threshold."
        ),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--factor-design", choices=("independent", "legacy_coupled"), default="independent")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    report = run_conformal_synthetic_gate(
        replicates=args.replicates, seed=args.seed, factor_design=args.factor_design
    )
    output = args.output or PROJECT_ROOT / "output" / (
        f"track_b_conformal_synthetic_gate_{args.factor_design}_{args.replicates}rep_seed{args.seed}.json"
    )
    if args.output is None and output.exists():
        raise FileExistsError(f"refusing to overwrite existing default output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report["aggregate"], indent=2, allow_nan=False))
    print(f"saved {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
