"""Prespecified, synthetic-only Track B forecast and uncertainty experiment.

The grid varies generating anisotropy and assumed boundary-noise amplitude. It keeps calibration cases
separate from evaluation cases. It reads no patient data and makes no clinical claims.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from src import fk_mask_inverse as fi
from src import fk_synthetic_validation as fv


SHAPE = (16, 16, 16)
H_MM = 2.0
RHO_PER_DAY = 0.04
D_ISOTROPIC_MM2_PER_DAY = 0.05
D_ANISOTROPIC_DIAGONAL_MM2_PER_DAY = (0.08, 0.05, 0.02)
TRAIN_DAYS = (0.0, 20.0, 40.0, 60.0)
TARGET_DAY = 80.0
CORRELATION_MM = 2.0
NOISE_LEVELS_MM = (0.0, 0.8)
TRUTHS = ("isotropic", "anisotropic")
ERROR_SCORE_NAMES = ("ensemble_disagreement", "training_spatial_residual",
                     "training_volume_residual_score", "composite_uncertainty")
HYBRID_ARM = "pde_shape_volume_trend"
FORECAST_ARMS = ("pde", "persistence", "volume_trend", HYBRID_ARM)


def build_split_manifest(replicates_per_split: int = 8, seed: int = 20261024) -> list[dict]:
    """Build balanced, disjoint calibration/evaluation case IDs before any forecast is run."""
    if isinstance(replicates_per_split, (bool, np.bool_)) or int(replicates_per_split) != replicates_per_split:
        raise ValueError("replicates_per_split must be a positive integer")
    replicates_per_split = int(replicates_per_split)
    if replicates_per_split < 1:
        raise ValueError("replicates_per_split must be a positive integer")
    if isinstance(seed, (bool, np.bool_)) or int(seed) != seed:
        raise ValueError("seed must be an integer")
    rows = []
    for condition_index, truth in enumerate(TRUTHS):
        for noise_index, sigma_mm in enumerate(NOISE_LEVELS_MM):
            for split_index, split in enumerate(("calibration", "evaluation")):
                for replicate in range(replicates_per_split):
                    case_seed = int(np.random.SeedSequence(
                        [int(seed), condition_index, noise_index, split_index, replicate]
                    ).generate_state(1)[0])
                    rows.append({
                        "case_id": f"{truth}_noise{sigma_mm:g}_{split}_{replicate:03d}",
                        "truth": truth,
                        "sigma_mm": float(sigma_mm),
                        "split": split,
                        "replicate": int(replicate),
                        "seed": case_seed,
                    })
    return rows


def _ball(shape, radius_voxels: float) -> np.ndarray:
    coords = np.ogrid[tuple(slice(0, n) for n in shape)]
    center = tuple((n - 1) / 2.0 for n in shape)
    return sum((axis - center_i) ** 2 for axis, center_i in zip(coords, center)) <= radius_voxels**2


def _truth_series(truth: str, initial: np.ndarray, domain: np.ndarray) -> list[np.ndarray]:
    days = [*TRAIN_DAYS, TARGET_DAY]
    if truth == "isotropic":
        future = fi.forward_masks(initial, domain, RHO_PER_DAY, D_ISOTROPIC_MM2_PER_DAY, days, H_MM)
    elif truth == "anisotropic":
        tensor = np.zeros(initial.shape + (6,), dtype=np.float32)
        tensor[..., :3] = D_ANISOTROPIC_DIAGONAL_MM2_PER_DAY
        future = fi.forward_masks_tensor(initial, domain, RHO_PER_DAY, tensor, days, H_MM)
    else:
        raise ValueError(f"unknown synthetic truth: {truth}")
    return [initial, *future]


def run_case(manifest_row: dict, *, bootstrap_draws: int = 4) -> dict:
    """Generate one case, fit only noisy training scans, and score one noisy future scan."""
    if manifest_row["truth"] not in TRUTHS or manifest_row["split"] not in ("calibration", "evaluation"):
        raise ValueError("manifest row has invalid truth or split")
    seed_sequence = np.random.SeedSequence(int(manifest_row["seed"]))
    children = seed_sequence.spawn(len(TRAIN_DAYS) + 1 + 1)
    rng = np.random.default_rng(children[0])
    radius = float(rng.choice((2.0, 2.5, 3.0, 3.5)))
    initial = _ball(SHAPE, radius)
    domain = np.ones(SHAPE, dtype=bool)
    truth_masks = _truth_series(manifest_row["truth"], initial, domain)
    sigma_mm = float(manifest_row["sigma_mm"])

    training = []
    for index, mask in enumerate(truth_masks[:len(TRAIN_DAYS)]):
        noisy = fv.perturb_boundary(
            mask, sigma_mm=sigma_mm, correlation_mm=CORRELATION_MM, spacing_mm=H_MM,
            seed=int(children[index + 1].generate_state(1)[0]),
        ) & domain
        if not noisy.any():
            raise ValueError(f"boundary noise emptied training mask for {manifest_row['case_id']}")
        training.append(noisy)
    target = fv.perturb_boundary(
        truth_masks[-1], sigma_mm=sigma_mm, correlation_mm=CORRELATION_MM, spacing_mm=H_MM,
        seed=int(children[-1].generate_state(1)[0]),
    ) & domain
    if not target.any():
        raise ValueError(f"boundary noise emptied target mask for {manifest_row['case_id']}")

    comparison = fv.compare_next_scan(
        training, TRAIN_DAYS, TARGET_DAY, target, h=H_MM, domain=domain
    )
    ensemble = fv.bootstrap_forecast_ensemble(
        training, TRAIN_DAYS, TARGET_DAY, h=H_MM, sigma_mm=sigma_mm,
        correlation_mm=CORRELATION_MM, n_draws=bootstrap_draws,
        seed=int(manifest_row["seed"]), domain=domain,
    )
    metrics = comparison["metrics"]
    pde_error = 1.0 - float(metrics["pde"]["dice"])
    ensemble_score = float(ensemble["disagreement"])
    spatial_score = float(comparison["training_spatial_residual"])
    volume_score = float(comparison["training_volume_residual_score"])
    return {
        **manifest_row,
        "radius_voxels": radius,
        "bootstrap_draws": int(bootstrap_draws),
        "pde_dice": float(metrics["pde"]["dice"]),
        "pde_error": pde_error,
        "pde_hd95_mm": float(metrics["pde"]["hd95_mm"]),
        "pde_symmetric_volume_error_pct": float(metrics["pde"]["spe_pct"]),
        "persistence_dice": float(metrics["persistence"]["dice"]),
        "persistence_error": 1.0 - float(metrics["persistence"]["dice"]),
        "persistence_hd95_mm": float(metrics["persistence"]["hd95_mm"]),
        "persistence_symmetric_volume_error_pct": float(metrics["persistence"]["spe_pct"]),
        "volume_trend_dice": float(metrics["volume_trend"]["dice"]),
        "volume_trend_error": 1.0 - float(metrics["volume_trend"]["dice"]),
        "volume_trend_hd95_mm": float(metrics["volume_trend"]["hd95_mm"]),
        "volume_trend_symmetric_volume_error_pct": float(metrics["volume_trend"]["spe_pct"]),
        f"{HYBRID_ARM}_dice": float(metrics[HYBRID_ARM]["dice"]),
        "hybrid_error": 1.0 - float(metrics[HYBRID_ARM]["dice"]),
        f"{HYBRID_ARM}_error": 1.0 - float(metrics[HYBRID_ARM]["dice"]),
        f"{HYBRID_ARM}_hd95_mm": float(metrics[HYBRID_ARM]["hd95_mm"]),
        f"{HYBRID_ARM}_symmetric_volume_error_pct": float(metrics[HYBRID_ARM]["spe_pct"]),
        "ensemble_disagreement": ensemble_score,
        "training_spatial_residual": spatial_score,
        "training_volume_residual_score": volume_score,
        "composite_uncertainty": float(np.mean([ensemble_score, spatial_score, volume_score])),
        "fit_rho": float(comparison["fit"]["rho"]),
        "fit_D_mm2_per_day": float(comparison["fit"]["D"]),
        "fit_rho_at_bound": bool(comparison["fit"]["rho_at_bound"]),
        "fit_D_at_bound": bool(comparison["fit"]["D_at_bound"]),
    }


def _mean_finite(values) -> float | None:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    return float(finite.mean()) if len(finite) else None


def summarize_rows(rows: list[dict], *, alpha: float = 0.2) -> dict:
    """Compare forecasts and calibrate within each generating condition.

    The aggregate coverage is descriptive. The finite-sample split-conformal statement applies within a
    condition only when calibration and future cases in that condition are exchangeable.
    """
    calibration = [row for row in rows if row["split"] == "calibration"]
    evaluation = [row for row in rows if row["split"] == "evaluation"]
    if not calibration or not evaluation:
        raise ValueError("need separate calibration and evaluation rows")
    if {row["case_id"] for row in calibration} & {row["case_id"] for row in evaluation}:
        raise ValueError("calibration and evaluation case IDs must be disjoint")

    def condition_key(row: dict) -> str:
        truth = row.get("truth")
        sigma_mm = row.get("sigma_mm")
        if truth not in TRUTHS or not isinstance(sigma_mm, (int, float)) or not math.isfinite(sigma_mm):
            raise ValueError("each row needs a known truth and finite sigma_mm condition")
        sigma_mm = float(sigma_mm)
        if sigma_mm not in NOISE_LEVELS_MM:
            raise ValueError(f"unsupported noise condition: {sigma_mm}")
        return f"{truth}|sigma_mm={sigma_mm:g}"

    calibration_by_condition = {}
    evaluation_by_condition = {}
    for destination, subset in ((calibration_by_condition, calibration), (evaluation_by_condition, evaluation)):
        for row in subset:
            destination.setdefault(condition_key(row), []).append(row)
    if set(calibration_by_condition) != set(evaluation_by_condition):
        raise ValueError("every generating condition needs both calibration and evaluation cases")

    forecast_comparison = {}
    for arm in FORECAST_ARMS:
        forecast_comparison[arm] = {
            "n_evaluation": len(evaluation),
            "mean_dice": _mean_finite([row[f"{arm}_dice"] for row in evaluation]),
            "mean_hd95_mm": _mean_finite([row[f"{arm}_hd95_mm"] for row in evaluation]),
            "mean_symmetric_volume_error_pct": _mean_finite(
                [row[f"{arm}_symmetric_volume_error_pct"] for row in evaluation]
            ),
        }

    uncertainty_by_condition = {}
    aggregate_predictions = {name: {"bounds": [], "errors": [], "n_calibration": 0}
                             for name in ERROR_SCORE_NAMES}
    for name in ERROR_SCORE_NAMES:
        aggregate_predictions[name]["n_calibration"] = len(calibration)
    for condition in sorted(calibration_by_condition):
        condition_calibration = calibration_by_condition[condition]
        condition_evaluation = evaluation_by_condition[condition]
        uncertainty_by_condition[condition] = {}
        for name in ERROR_SCORE_NAMES:
            calibrated = fv.conformal_error_upper_bounds(
                [row[name] for row in condition_calibration],
                [row["hybrid_error"] for row in condition_calibration],
                [row[name] for row in condition_evaluation],
                alpha=alpha,
            )
            bounds = calibrated["bounds"]
            errors = np.asarray([row["hybrid_error"] for row in condition_evaluation], dtype=float)
            uncertainty_by_condition[condition][name] = {
                "n_calibration": len(condition_calibration),
                "n_evaluation": len(condition_evaluation),
                "calibration_quantile": calibrated["quantile"],
                "evaluation_coverage": float(np.mean(errors <= bounds)),
                "nominal_coverage": float(1.0 - alpha),
                "mean_upper_bound": float(np.mean(bounds)),
                "mean_hybrid_error": float(errors.mean()),
                "target_forecast": HYBRID_ARM,
            }
            aggregate_predictions[name]["bounds"].extend(bounds.tolist())
            aggregate_predictions[name]["errors"].extend(errors.tolist())

    uncertainty = {}
    for name, aggregate in aggregate_predictions.items():
        bounds = np.asarray(aggregate["bounds"], dtype=float)
        errors = np.asarray(aggregate["errors"], dtype=float)
        uncertainty[name] = {
            "n_calibration": aggregate["n_calibration"],
            "n_evaluation": int(len(errors)),
            "evaluation_coverage": float(np.mean(errors <= bounds)),
            "nominal_coverage": float(1.0 - alpha),
            "mean_upper_bound": float(np.mean(bounds)),
            "mean_hybrid_error": float(errors.mean()),
            "target_forecast": HYBRID_ARM,
            "coverage_interpretation": "descriptive aggregate of condition-specific bounds",
        }
    return {"alpha": float(alpha), "forecast_comparison": forecast_comparison,
            "uncertainty": uncertainty, "uncertainty_by_condition": uncertainty_by_condition,
            "pde_minus_comparator_dice_error": {
                arm: _mean_finite([row["pde_error"] - row[f"{arm}_error"] for row in evaluation])
                for arm in ("persistence", "volume_trend", HYBRID_ARM)
            }}


def _json_safe(value):
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def run_grid(*, replicates_per_split: int = 8, bootstrap_draws: int = 4,
             seed: int = 20261024, alpha: float = 0.2) -> tuple[list[dict], dict]:
    """Run the fixed 2 x 2 synthetic grid, with equal calibration/evaluation counts per cell."""
    if isinstance(bootstrap_draws, (bool, np.bool_)) or int(bootstrap_draws) != bootstrap_draws or bootstrap_draws < 2:
        raise ValueError("bootstrap_draws must be an integer of at least two")
    manifest = build_split_manifest(replicates_per_split=replicates_per_split, seed=seed)
    rows = [run_case(row, bootstrap_draws=int(bootstrap_draws)) for row in manifest]
    summary = summarize_rows(rows, alpha=alpha)
    summary["design"] = {
        "truths": list(TRUTHS),
        "noise_levels_mm": list(NOISE_LEVELS_MM),
        "correlation_mm": CORRELATION_MM,
        "shape_voxels": list(SHAPE),
        "spacing_mm": H_MM,
        "rho_per_day": RHO_PER_DAY,
        "hybrid_method": "PDE spatial forecast projected to the independent exponential volume-trend size",
        "D_isotropic_mm2_per_day": D_ISOTROPIC_MM2_PER_DAY,
        "D_anisotropic_diagonal_mm2_per_day": list(D_ANISOTROPIC_DIAGONAL_MM2_PER_DAY),
        "training_days": list(TRAIN_DAYS),
        "target_day": TARGET_DAY,
        "replicates_per_split_per_cell": int(replicates_per_split),
        "bootstrap_draws_per_case": int(bootstrap_draws),
        "seed": int(seed),
        "cases_total": len(rows),
        "case_status": "synthetic only; no patient data",
        "noise_status": "assumed spatially correlated boundary perturbation; not empirically measured",
    }
    return rows, summary


def write_outputs(rows: list[dict], summary: dict, output_dir: Path) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path, json_path = output_dir / "case_metrics.csv", output_dir / "summary.json"
    if csv_path.exists() or json_path.exists():
        raise FileExistsError(f"refusing to overwrite existing scenario-grid outputs in {output_dir}")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(_json_safe(summary), indent=2, allow_nan=False) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--replicates-per-split", type=int, default=8)
    parser.add_argument("--bootstrap-draws", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20261024)
    parser.add_argument("--alpha", type=float, default=0.2)
    args = parser.parse_args(argv)
    rows, summary = run_grid(replicates_per_split=args.replicates_per_split,
                             bootstrap_draws=args.bootstrap_draws, seed=args.seed, alpha=args.alpha)
    write_outputs(rows, summary, args.output_dir)
    print(json.dumps(_json_safe(summary), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
