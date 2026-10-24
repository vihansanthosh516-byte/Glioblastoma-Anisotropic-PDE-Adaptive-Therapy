"""Outer-patient evaluation of the Track B hybrid on a controlled synthetic cohort.

All inputs are generated here. The evaluator uses no patient data. The simulated segmentation-error
covariance is deliberately supplied at its known generating value, so this is a favorable method check and
does not validate noise estimation or performance on MRI cohorts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from src import fk_mask_inverse as fi
from src.fk_boundary_calibration import run_cross_fitted_boundary_calibration
from src.fk_grid_likelihood import volume_match_concentration
from src.fk_hybrid_forecast import estimate_log_volume_slope
from src.fk_patient_cv import PatientSeries, run_outer_patient_cv
from src.fk_spatial_uncertainty import (
    correlated_boundary_perturbation,
    energy_score_mask,
    fixed_boundary_score_pairs,
    marginal_brier_score,
    spatial_variogram_score,
    volume_crps,
)
from src.fk_synthetic_validation import ensemble_boundary_disagreement
from src.fk_volume_model import fit_log_volume_variance, select_shrinkage_time_constant, shrink_log_growth_slope


PROJECT_ROOT = Path(__file__).resolve().parent.parent
SHAPE = (16, 16, 16)
SPACING_MM = 2.0
TRAINING_DAYS = (0.0, 15.0, 35.0)
DIFFUSION_TRUTHS = ("isotropic", "anisotropic")
VOLUME_TRENDS = ("growth", "shrinkage")
BOUNDARY_TRUTH = {"sigma_b_mm": 0.6, "tau_b_mm": 0.2, "ell_mm": 4.0}
# Fixed working value used only to construct cross-fitted forecasts before residual covariance is estimated.
# It is deliberately distinct from the synthetic generator's covariance and is not a claimed scanner-noise value.
BOUNDARY_CALIBRATION_START = {"sigma_b_mm": 0.8, "tau_b_mm": 0.3, "ell_mm": 6.0}


def build_synthetic_manifest(replicates_per_condition: int = 25, n_folds: int = 5,
                             seed: int = 20261024) -> list[dict]:
    """Build a balanced four-condition manifest with patient-level outer folds."""
    if (isinstance(replicates_per_condition, (bool, np.bool_))
            or int(replicates_per_condition) != replicates_per_condition
            or replicates_per_condition < 2):
        raise ValueError("replicates_per_condition must be an integer of at least two")
    if isinstance(n_folds, (bool, np.bool_)) or int(n_folds) != n_folds or n_folds < 2:
        raise ValueError("n_folds must be an integer of at least two")
    if int(replicates_per_condition) % int(n_folds):
        raise ValueError("replicates_per_condition must be divisible by n_folds for balanced splits")
    if isinstance(seed, (bool, np.bool_)) or int(seed) != seed:
        raise ValueError("seed must be an integer")
    replicates_per_condition, n_folds = int(replicates_per_condition), int(n_folds)
    rows = []
    condition_index = 0
    for diffusion_truth in DIFFUSION_TRUTHS:
        for volume_trend in VOLUME_TRENDS:
            for replicate in range(replicates_per_condition):
                case_seed = int(np.random.SeedSequence(
                    [int(seed), condition_index, replicate]
                ).generate_state(1)[0])
                rows.append({
                    "case_id": f"syn-{diffusion_truth}-{volume_trend}-{replicate:03d}",
                    "diffusion_truth": diffusion_truth,
                    "volume_trend": volume_trend,
                    "replicate": replicate,
                    "fold": replicate % n_folds,
                    "seed": case_seed,
                })
            condition_index += 1
    return rows


def _ball(radius_voxels: float, shape=SHAPE) -> np.ndarray:
    coords = np.ogrid[tuple(slice(0, n) for n in shape)]
    center = tuple((n - 1) / 2.0 for n in shape)
    return sum((axis - center_i) ** 2 for axis, center_i in zip(coords, center)) <= radius_voxels**2


def generate_synthetic_series(manifest_row: dict) -> tuple[PatientSeries, dict]:
    """Generate four 3D visits with PDE-driven shape, separate growth/shrinkage volume, and correlated error."""
    diffusion_truth = manifest_row.get("diffusion_truth")
    volume_trend = manifest_row.get("volume_trend")
    if diffusion_truth not in DIFFUSION_TRUTHS or volume_trend not in VOLUME_TRENDS:
        raise ValueError("manifest row has an unsupported generating condition")
    seed = manifest_row.get("seed")
    if isinstance(seed, (bool, np.bool_)) or int(seed) != seed:
        raise ValueError("manifest row seed must be an integer")
    rng = np.random.default_rng(int(seed))
    replicate = int(manifest_row["replicate"])
    interval_cycle = (15.0, 20.0, 25.0)
    gaps = (interval_cycle[replicate % 3], interval_cycle[(replicate + 1) % 3],
            interval_cycle[(replicate + 2) % 3])
    days = np.asarray([0.0, gaps[0], gaps[0] + gaps[1], sum(gaps)], dtype=float)
    domain = np.ones(SHAPE, dtype=bool)
    initial = _ball(float(rng.choice((3.0, 3.5, 4.0))))
    tensor = np.zeros(SHAPE + (6,), dtype=np.float32)
    if diffusion_truth == "isotropic":
        tensor[..., :3] = 0.05
    else:
        high_axis = replicate % 3
        diagonal = np.asarray((0.08, 0.035, 0.02), dtype=np.float32)
        diagonal = np.roll(diagonal, high_axis)
        tensor[..., :3] = diagonal

    concentrations = [initial.astype(np.float32), *fi.forward_concentrations_tensor(
        initial, domain, 0.04, tensor, days, SPACING_MM
    )]
    rate_magnitude = float(rng.uniform(0.004, 0.010))
    sign = 1.0 if volume_trend == "growth" else -1.0
    slope = sign * rate_magnitude
    initial_volume = int(initial.sum())
    masks = []
    for scan_index, (day, concentration) in enumerate(zip(days, concentrations)):
        requested_volume = int(np.rint(initial_volume * np.exp(slope * day)))
        requested_volume = int(np.clip(requested_volume, 8, int(domain.sum()) - 1))
        clean_mask = volume_match_concentration(
            concentration, domain, target_voxels=requested_volume
        )
        observed = correlated_boundary_perturbation(
            clean_mask, **BOUNDARY_TRUTH, spacing_mm=SPACING_MM,
            seed=int(np.random.SeedSequence([int(seed), scan_index, 71]).generate_state(1)[0]),
            n_anchors=24,
        ) & domain
        if not observed.any() or observed.all():
            raise ValueError(f"synthetic boundary error made {manifest_row['case_id']} unusable")
        masks.append(observed)

    record = PatientSeries(str(manifest_row["case_id"]), days, tuple(masks))
    truth = {
        "diffusion_truth": diffusion_truth,
        "volume_trend": volume_trend,
        "true_log_volume_slope_per_day": slope,
        "true_tensor_diagonal_mm2_per_day": tuple(float(value) for value in tensor[0, 0, 0, :3]),
        "boundary_covariance_truth": dict(BOUNDARY_TRUTH),
    }
    return record, truth


def fit_training_nuisance(training_records) -> dict:
    """Fit volume-population slope, shrinkage, and variance using only the supplied training patients.

    Boundary covariance and the parameter prior are fixed to the known synthetic generator values. The
    covariance is an oracle input here; real use still needs outer-training cross-fitted residual estimation.
    """
    records = tuple(training_records)
    if len(records) < 3:
        raise ValueError("need at least three outer-training patients")
    ids = [record.patient_id for record in records]
    if len(set(ids)) != len(ids):
        raise ValueError("training records must have unique patient IDs")
    patient_slopes = {}
    for record in records:
        volumes = np.asarray([np.asarray(mask, dtype=bool).sum() for mask in record.masks], dtype=float)
        patient_slopes[record.patient_id] = estimate_log_volume_slope(volumes, record.days)

    validation_cases = []
    for record in records:
        other_slopes = [slope for patient_id, slope in patient_slopes.items()
                        if patient_id != record.patient_id]
        if not other_slopes:
            raise ValueError("leave-one-patient-out validation needs other training patients")
        leave_one_out_population = float(np.mean(other_slopes))
        volumes = np.asarray([np.asarray(mask, dtype=bool).sum() for mask in record.masks], dtype=float)
        for target_index in range(2, len(record.days)):
            history_days = np.asarray(record.days[:target_index], dtype=float)
            history_volumes = volumes[:target_index]
            validation_cases.append({
                "patient_slope": estimate_log_volume_slope(history_volumes, history_days),
                "population_slope": leave_one_out_population,
                "observation_gap_days": float(history_days[-1] - history_days[0]),
                "horizon_days": float(record.days[target_index] - history_days[-1]),
                "last_volume": float(history_volumes[-1]),
                "target_volume": float(volumes[target_index]),
            })
    selected_tau = select_shrinkage_time_constant(validation_cases)

    validation_horizons, validation_residuals = [], []
    for case in validation_cases:
        shrunk = shrink_log_growth_slope(
            case["patient_slope"], case["population_slope"],
            observation_gap_days=case["observation_gap_days"],
            time_constant_days=selected_tau["time_constant_days"],
        )
        prediction = np.log(case["last_volume"]) + shrunk * case["horizon_days"]
        validation_horizons.append(case["horizon_days"])
        validation_residuals.append(np.log(case["target_volume"]) - prediction)
    variance_fit = fit_log_volume_variance(validation_horizons, validation_residuals)
    return {
        "training_patient_ids": tuple(sorted(ids)),
        "population_log_volume_slope_per_day": float(np.mean(list(patient_slopes.values()))),
        "shrinkage_time_constant_days": float(selected_tau["time_constant_days"]),
        "shrinkage_selection_errors": selected_tau["candidate_mse_log_volume"],
        "volume_variance_fit": variance_fit,
        "prior": {"mean": np.log([0.04, 0.05]), "covariance": np.eye(2) * 0.25},
        "boundary_covariance": dict(BOUNDARY_TRUTH),
        "boundary_covariance_source": "known simulation generator (oracle)",
    }


def _seed_for_patient(patient_id: str) -> int:
    return int.from_bytes(hashlib.sha256(patient_id.encode("utf-8")).digest()[:4], "little")


def _forecast_callback(history, nuisance, *, forecast_function, forecast_config):
    shape = history.masks[-1].shape
    result = forecast_function(
        history.masks, history.days, history.forecast_day,
        prior=nuisance["prior"], h=SPACING_MM,
        boundary_covariance=nuisance["boundary_covariance"],
        population_log_volume_slope_per_day=nuisance["population_log_volume_slope_per_day"],
        shrinkage_time_constant_days=nuisance["shrinkage_time_constant_days"],
        volume_variance_fit=nuisance["volume_variance_fit"],
        domain=np.ones(shape, dtype=bool),
        rho_bounds=(0.01, 0.1), D_bounds=(0.01, 0.1),
        n_coarse=forecast_config["n_coarse"], n_refine=forecast_config["n_refine"],
        n_draws=forecast_config["n_draws"], n_anchors=forecast_config["n_anchors"],
        seed=_seed_for_patient(history.patient_id),
    )
    return {"forecast": result, "history_mask": history.masks[-1]}


def _score_callback(prediction, target) -> dict:
    forecast = prediction["forecast"]
    hybrid_draws = np.asarray(forecast["mask_draws"], dtype=bool)
    baseline_draws = np.asarray(forecast["no_pde_ranker_mask_draws"], dtype=bool)
    truth = np.asarray(target.mask, dtype=bool)
    history_mask = np.asarray(prediction["history_mask"], dtype=bool)
    if (hybrid_draws.ndim != 4 or baseline_draws.shape != hybrid_draws.shape
            or truth.shape != hybrid_draws.shape[1:] or history_mask.shape != truth.shape):
        raise ValueError("forecast draws, held-out mask, and history mask must have matching 3D shapes")
    design = fixed_boundary_score_pairs(
        history_mask=history_mask, spacing_mm=SPACING_MM, n_voxels=128,
        band_mm=10.0, max_pair_distance_mm=10.0,
    )
    hybrid_consensus = hybrid_draws.mean(axis=0) >= 0.5
    baseline_consensus = baseline_draws.mean(axis=0) >= 0.5
    hybrid_energy = energy_score_mask(hybrid_draws, truth, design["voxel_indices"])
    baseline_energy = energy_score_mask(baseline_draws, truth, design["voxel_indices"])
    hybrid_variogram = spatial_variogram_score(hybrid_draws, truth, design)["score"]
    baseline_variogram = spatial_variogram_score(baseline_draws, truth, design)["score"]
    hybrid_brier = marginal_brier_score(hybrid_draws, truth, design["voxel_indices"])
    baseline_brier = marginal_brier_score(baseline_draws, truth, design["voxel_indices"])
    hybrid_volume_crps = volume_crps(hybrid_draws.reshape(len(hybrid_draws), -1).sum(axis=1), truth.sum())
    baseline_volume_crps = volume_crps(baseline_draws.reshape(len(baseline_draws), -1).sum(axis=1), truth.sum())
    from src.forecast_metrics import all_metrics
    hybrid_dice = float(all_metrics(hybrid_consensus, truth, spacing=SPACING_MM)["dice"])
    baseline_dice = float(all_metrics(baseline_consensus, truth, spacing=SPACING_MM)["dice"])
    return {
        "hybrid_energy_score": float(hybrid_energy),
        "no_pde_energy_score": float(baseline_energy),
        "energy_score_difference_hybrid_minus_no_pde": float(hybrid_energy - baseline_energy),
        "hybrid_variogram_score": float(hybrid_variogram),
        "no_pde_variogram_score": float(baseline_variogram),
        "hybrid_brier_score": float(hybrid_brier),
        "no_pde_brier_score": float(baseline_brier),
        "hybrid_volume_crps_voxels": float(hybrid_volume_crps),
        "no_pde_volume_crps_voxels": float(baseline_volume_crps),
        "hybrid_consensus_dice": hybrid_dice,
        "no_pde_consensus_dice": baseline_dice,
        "hybrid_minus_no_pde_dice": float(hybrid_dice - baseline_dice),
        "hybrid_lower_energy": bool(hybrid_energy < baseline_energy),
        "hybrid_boundary_disagreement": float(ensemble_boundary_disagreement(hybrid_draws)),
        "target_volume_voxels": int(truth.sum()),
        "forecast_day": float(target.day),
    }


def _summarize(rows: list[dict]) -> dict:
    def summary(subset):
        differences = np.asarray(
            [row["scores"]["energy_score_difference_hybrid_minus_no_pde"] for row in subset], dtype=float
        )
        return {
            "n_patients": int(len(subset)),
            "mean_energy_score_difference_hybrid_minus_no_pde": float(differences.mean()),
            "median_energy_score_difference_hybrid_minus_no_pde": float(np.median(differences)),
            "hybrid_lower_energy_count": int(np.sum(differences < 0)),
            "mean_hybrid_consensus_dice": float(np.mean([r["scores"]["hybrid_consensus_dice"] for r in subset])),
            "mean_no_pde_consensus_dice": float(np.mean([r["scores"]["no_pde_consensus_dice"] for r in subset])),
            "mean_hybrid_variogram_score": float(np.mean([r["scores"]["hybrid_variogram_score"] for r in subset])),
            "mean_no_pde_variogram_score": float(np.mean([r["scores"]["no_pde_variogram_score"] for r in subset])),
            "mean_hybrid_volume_crps_voxels": float(np.mean([r["scores"]["hybrid_volume_crps_voxels"] for r in subset])),
            "mean_no_pde_volume_crps_voxels": float(np.mean([r["scores"]["no_pde_volume_crps_voxels"] for r in subset])),
        }

    by_condition = {}
    for row in rows:
        key = f"{row['scores']['diffusion_truth']}|{row['scores']['volume_trend']}"
        by_condition.setdefault(key, []).append(row)
    return {
        "n_patients": len(rows),
        "overall": summary(rows),
        "by_condition": {key: summary(value) for key, value in sorted(by_condition.items())},
    }


def run_synthetic_gate(*, replicates_per_condition: int = 25, n_folds: int = 5,
                       seed: int = 20261024, forecast_function=None,
                       forecast_config: dict | None = None,
                       forecast_nuisance_overrides: dict | None = None,
                       boundary_calibration_mode: str = "cross_fitted",
                       inner_n_folds: int = 2) -> dict:
    """Run the hybrid on synthetic, patient-held-out histories and score matched-size ablations.

    By default, each outer fold estimates boundary covariance from inner patient-held-out residuals drawn
    only from that fold's training patients. ``oracle`` mode is retained only to reproduce the favorable
    historical simulation; it supplies the generator's known covariance.
    """
    from src import fk_hybrid_forecast as hybrid
    forecast_function = forecast_function or hybrid.hybrid_posterior_predictive_forecast
    forecast_config = {"n_coarse": 5, "n_refine": 3, "n_draws": 16, "n_anchors": 24,
                       **(forecast_config or {})}
    if boundary_calibration_mode not in {"cross_fitted", "oracle"}:
        raise ValueError("boundary_calibration_mode must be 'cross_fitted' or 'oracle'")
    if (isinstance(inner_n_folds, (bool, np.bool_)) or int(inner_n_folds) != inner_n_folds
            or inner_n_folds < 2):
        raise ValueError("inner_n_folds must be an integer of at least two")
    inner_n_folds = int(inner_n_folds)
    overrides = {} if forecast_nuisance_overrides is None else dict(forecast_nuisance_overrides)
    allowed_overrides = {"boundary_covariance", "prior"}
    if set(overrides) - allowed_overrides:
        raise ValueError(f"unsupported forecast nuisance overrides: {sorted(set(overrides) - allowed_overrides)}")
    if boundary_calibration_mode == "cross_fitted" and "boundary_covariance" in overrides:
        raise ValueError("boundary covariance overrides require oracle mode; cross-fitted mode estimates it from training patients")
    model_covariance = dict(overrides.get("boundary_covariance", BOUNDARY_TRUTH))
    required_covariance = {"sigma_b_mm", "tau_b_mm", "ell_mm"}
    if set(model_covariance) != required_covariance:
        raise ValueError("boundary covariance override must specify sigma_b_mm, tau_b_mm, and ell_mm")
    try:
        model_covariance = {key: float(value) for key, value in model_covariance.items()}
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("boundary covariance override values must be numeric") from exc
    if (not all(np.isfinite(value) for value in model_covariance.values())
            or model_covariance["sigma_b_mm"] < 0 or model_covariance["tau_b_mm"] < 0
            or model_covariance["ell_mm"] <= 0):
        raise ValueError("boundary covariance override parameters are outside their valid ranges")
    model_prior = overrides.get("prior", {"mean": np.log([0.04, 0.05]), "covariance": np.eye(2) * 0.25})
    if not isinstance(model_prior, dict):
        raise ValueError("parameter prior override must be a mapping")
    prior_mean = np.asarray(model_prior.get("mean"), dtype=float)
    prior_covariance = np.asarray(model_prior.get("covariance"), dtype=float)
    if (prior_mean.shape != (2,) or prior_covariance.shape != (2, 2)
            or not np.isfinite(prior_mean).all() or not np.isfinite(prior_covariance).all()
            or not np.allclose(prior_covariance, prior_covariance.T)
            or np.linalg.eigvalsh(prior_covariance).min() <= 0):
        raise ValueError("parameter prior override must have a finite 2-vector mean and positive-definite covariance")
    model_prior = {"mean": prior_mean, "covariance": prior_covariance}

    boundary_calibration_fits = []

    def fit_outer_training_nuisance(train_records):
        fitted = fit_training_nuisance(train_records)
        fitted["prior"] = model_prior
        if boundary_calibration_mode == "oracle":
            fitted["boundary_covariance"] = dict(model_covariance)
            fitted["boundary_covariance_source"] = (
                "prespecified sensitivity override" if overrides
                else "known simulation generator (oracle)"
            )
            return fitted

        if len(train_records) - int(np.ceil(len(train_records) / inner_n_folds)) < 3:
            raise ValueError("each inner calibration training split must contain at least three patients")

        def fit_inner_nuisance(inner_training_records):
            inner = fit_training_nuisance(inner_training_records)
            inner["prior"] = model_prior
            # The covariance enters the PDE shape likelihood, so it must remain positive definite. This
            # prespecified working value is not the generator truth and is replaced by cross-fit residuals.
            inner["boundary_covariance"] = dict(BOUNDARY_CALIBRATION_START)
            return inner

        def forecast_for_calibration(history, inner_nuisance):
            prediction = _forecast_callback(
                history, inner_nuisance, forecast_function=forecast_function, forecast_config=forecast_config,
            )
            forecast = prediction["forecast"]
            if "pre_boundary_probability" not in forecast:
                raise ValueError("cross-fitted boundary calibration requires pre_boundary_probability")
            return {"pre_boundary_probability": forecast["pre_boundary_probability"]}

        training_ids = sorted(record.patient_id for record in train_records)
        calibration_seed = int.from_bytes(
            hashlib.sha256("\0".join(training_ids).encode("utf-8")).digest()[:4], "little"
        )
        calibration = run_cross_fitted_boundary_calibration(
            train_records, fit_inner_nuisance, forecast_for_calibration,
            n_folds=inner_n_folds, seed=calibration_seed,
            minimum_history=3, spacing_mm=SPACING_MM, n_anchors=forecast_config["n_anchors"],
        )
        covariance_fit = calibration["covariance_fit"]
        estimated_covariance = {
            key: float(covariance_fit[key]) for key in ("sigma_b_mm", "tau_b_mm", "ell_mm")
        }
        fitted["boundary_covariance"] = estimated_covariance
        fitted["boundary_covariance_source"] = "inner patient-held-out residuals from outer-training patients"
        boundary_calibration_fits.append({
            "n_outer_training_patients": len(train_records),
            "n_calibration_patients": int(covariance_fit["n_calibration_patients"]),
            "n_calibration_origins": int(covariance_fit["n_calibration_origins"]),
            "estimate": estimated_covariance,
            "selection_boundary": bool(covariance_fit["selection_boundary"]),
            "weighted_rmse_mm2": float(covariance_fit["candidate_fits"][covariance_fit["ell_mm"]]["weighted_rmse_mm2"]),
            "fit_scope": covariance_fit["scope"],
        })
        return fitted

    manifest = build_synthetic_manifest(replicates_per_condition, n_folds, seed)
    generated = [generate_synthetic_series(row) for row in manifest]
    records = [record for record, _ in generated]
    truth_by_patient = {record.patient_id: truth for record, truth in generated}
    folds = {row["case_id"]: int(row["fold"]) for row in manifest}
    started = time.perf_counter()
    result = run_outer_patient_cv(
        records, folds, fit_outer_training_nuisance,
        lambda history, nuisance: _forecast_callback(
            history, nuisance, forecast_function=forecast_function, forecast_config=forecast_config,
        ),
        _score_callback, n_folds=n_folds, minimum_history=3,
    )
    elapsed = time.perf_counter() - started
    for row in result["rows"]:
        truth = truth_by_patient[row["patient_id"]]
        row["scores"].update({
            "diffusion_truth": truth["diffusion_truth"],
            "volume_trend": truth["volume_trend"],
        })
    boundary_search_boundary_folds = sum(
        fit["selection_boundary"] for fit in boundary_calibration_fits
    )
    if boundary_calibration_mode == "oracle":
        boundary_calibration_status = "not_run_oracle"
    elif boundary_search_boundary_folds:
        boundary_calibration_status = "inconclusive_search_boundary"
    else:
        boundary_calibration_status = "all_folds_interior"
    return {
        "manifest": manifest,
        "rows": result["rows"],
        "fold_audit": result["fold_audit"],
        "summary": _summarize(result["rows"]),
        "metadata": {
            "scope": "synthetic-only outer-patient rolling-origin simulation",
            "n_patients": len(records),
            "conditions": [f"{d}|{v}" for d in DIFFUSION_TRUTHS for v in VOLUME_TRENDS],
            "solver_spacing_mm": SPACING_MM,
            "shape_grid": SHAPE,
            "boundary_covariance_truth": dict(BOUNDARY_TRUTH),
            "boundary_calibration_mode": boundary_calibration_mode,
            "boundary_covariance_model_input": dict(model_covariance) if boundary_calibration_mode == "oracle" else None,
            "boundary_covariance_is_oracle_input": boundary_calibration_mode == "oracle" and model_covariance == BOUNDARY_TRUTH,
            "boundary_covariance_crossfit_start_input": (dict(BOUNDARY_CALIBRATION_START)
                                                           if boundary_calibration_mode == "cross_fitted" else None),
            "boundary_covariance_search_boundary_fold_count": int(boundary_search_boundary_folds),
            "boundary_covariance_calibration_status": boundary_calibration_status,
            "boundary_covariance_status_rule": (
                "inconclusive_search_boundary if any outer-fold length-scale estimate is at either candidate-grid edge; "
                "all_folds_interior is diagnostic only and does not establish calibration"
            ),
            "boundary_covariance_estimates_by_outer_fold": boundary_calibration_fits,
            "parameter_prior_model_input": {"mean": prior_mean.tolist(), "covariance": prior_covariance.tolist()},
            "parameter_prior_is_fixed_across_cases": True,
            "sensitivity_override_keys": sorted(overrides),
            "volume_population_slope_tau_and_variance_fit_scope": "outer-training patients only; tau via leave-one-patient-out rolling origins",
            "primary_comparison": "paired spatial energy score, hybrid PDE shape vs no-PDE signed-distance ranker, same sampled volumes",
            "secondary_measures": ["Dice", "spatial variogram score", "marginal Brier", "volume CRPS"],
            "elapsed_seconds_forecasting_and_scoring": float(elapsed),
            "patient_images_read": False,
            "interpretation_limit": "controlled simulation; not empirical MRI evidence or calibrated scanner-noise validation",
        },
        "scope": result["scope"],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates-per-condition", type=int, default=25)
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20261024)
    parser.add_argument("--boundary-calibration-mode", choices=("cross_fitted", "oracle"), default="cross_fitted")
    parser.add_argument("--inner-n-folds", type=int, default=2)
    parser.add_argument("--output", type=Path,
                        default=PROJECT_ROOT / "output" / "track_b_full_hybrid_gate_crossfit_20261024.json")
    args = parser.parse_args(argv)
    result = run_synthetic_gate(
        replicates_per_condition=args.replicates_per_condition,
        n_folds=args.n_folds, seed=args.seed,
        boundary_calibration_mode=args.boundary_calibration_mode,
        inner_n_folds=args.inner_n_folds,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    print(f"saved {args.output}")
    print(json.dumps(result["summary"], indent=2))
    print(json.dumps(result["metadata"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
