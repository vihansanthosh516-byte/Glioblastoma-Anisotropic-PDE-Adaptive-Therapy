"""Gated patient-level MU development runner for the fixed day-90 PDE protocol.

This runner is an exploratory development workflow on previously examined MU
data. The caller must first obtain a ``MUDevelopmentCohort`` from the loader,
which refuses mask access until registration and target-support evidence passes.
Boundary calibration starts from a caller-prespecified covariance and is fit
inside each outer training fold; the caller must lock that seed before scoring.
"""
from __future__ import annotations

from collections.abc import Mapping
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
from types import MappingProxyType

import numpy as np
import scipy

from src.fk_boundary_calibration import run_cross_fitted_boundary_calibration
from src.fk_mu_development_data import MUDevelopmentCohort, is_loader_issued_cohort
from src.fk_hybrid_forecast import hybrid_posterior_predictive_forecast
from src.fk_mu_scoring import make_first_third_score_callback
from src.fk_mu_training_nuisance import fit_mu_training_volume_nuisance
from src.fk_patient_cv import ForecastFailure, PatientHistory, run_first_third_90day_patient_cv
from src.fk_precision_gate import (
    ExternalEligibilityEvidence,
    is_validated_external_eligibility,
    precision_gate_from_patient_differences,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_ROOT = PROJECT_ROOT / "output"
# Match the existing synthetic-gate working start; this is an initialization, not measured noise or fitted truth.
LOCKED_MU_INITIAL_BOUNDARY_COVARIANCE = MappingProxyType({
    "sigma_b_mm": 0.8,
    "tau_b_mm": 0.3,
    "ell_mm": 6.0,
})


_EXPECTED_MODEL_FAILURES = {
    "forecast volume reaches computational-domain capacity; enlarge the domain",
    "volume-matched forecast has no usable boundary",
}


def _stable_seed(base_seed: int, patient_id: str) -> int:
    digest = hashlib.sha256(f"{int(base_seed)}\0{patient_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "little", signed=False)


def write_mu_development_result(result: Mapping, output_path: str | Path) -> dict:
    """Write one no-overwrite JSON result with code, protocol, data, and environment hashes."""
    if (not isinstance(result, Mapping) or type(result.get("schema_version")) is not int
            or result.get("schema_version") != 1):
        raise ValueError("result must be a schema-1 MU development result")
    path = Path(output_path).resolve()
    output_root = OUTPUT_ROOT.resolve()
    if path != output_root and output_root not in path.parents:
        raise PermissionError("MU development results must be written inside the project output directory")
    if path.suffix.lower() != ".json":
        raise ValueError("MU development result path must end in .json")

    code_files = (
        "src/fk_mu_patient_runner.py",
        "src/fk_mu_development_data.py",
        "src/fk_mu_training_nuisance.py",
        "src/fk_mu_scoring.py",
        "src/fk_hybrid_forecast.py",
        "src/fk_grid_likelihood.py",
        "src/fk_joint_posterior.py",
        "src/fk_mask_inverse.py",
        "src/fk_volume_model.py",
        "src/fk_spatial_uncertainty.py",
        "src/fk_precision_gate.py",
        "src/fk_boundary_calibration.py",
        "src/fk_patient_cv.py",
        "run_improved_aniso.py",
    )
    code_hashes = {}
    for relative_path in code_files:
        code_path = PROJECT_ROOT / relative_path
        code_hashes[relative_path] = hashlib.sha256(code_path.read_bytes()).hexdigest()
    protocol_path = PROJECT_ROOT / "paper/track_b_pde_external_protocol_draft_2026-10-28.md"
    protocol_hash = hashlib.sha256(protocol_path.read_bytes()).hexdigest()
    nested_runs = result.get("runs")
    if isinstance(nested_runs, Mapping):
        expected_run_names = {"primary", "one_decade_narrower", "one_decade_wider"}
        if set(nested_runs) != expected_run_names or any(
                not isinstance(run, Mapping) for run in nested_runs.values()):
            raise ValueError("protocol-complete result must contain all three locked prior-support runs")
        boundary_gap_run = result.get("boundary_gap_sensitivity_run")
        if not isinstance(boundary_gap_run, Mapping):
            raise ValueError("protocol-complete result must contain the central-gap sensitivity run")
        primary_run = nested_runs.get("primary")
        if not isinstance(primary_run, Mapping):
            raise ValueError("protocol-complete result must contain a primary run")
        cohort_audit = primary_run.get("cohort_audit", {})
        run_settings = {
            "runs": {
                name: {
                    "protocol": run.get("protocol", {}),
                    "random_seeds": run.get("random_seeds", {}),
                }
                for name, run in nested_runs.items()
            },
            "boundary_gap_sensitivity_run": {
                "protocol": boundary_gap_run.get("protocol", {}),
                "random_seeds": boundary_gap_run.get("random_seeds", {}),
            },
        }
    else:
        cohort_audit = result.get("cohort_audit", {})
        run_settings = {
            "protocol": result.get("protocol", {}),
            "random_seeds": result.get("random_seeds", {}),
        }
    data_hashes = cohort_audit.get("source_hashes", {}) if isinstance(cohort_audit, Mapping) else {}
    if (not isinstance(data_hashes, Mapping) or not data_hashes
            or any(not isinstance(name, str) or not name.strip()
                   or not isinstance(digest, str) or len(digest) != 64
                   or any(char not in "0123456789abcdef" for char in digest.lower())
                   for name, digest in data_hashes.items())):
        raise ValueError("MU result cohort audit must contain nonempty SHA-256 source hashes")
    environment = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    run_spec = {
        "run_settings": run_settings,
        "source_hashes": dict(data_hashes),
        "code_sha256": code_hashes,
        "protocol_sha256": protocol_hash,
        "environment": environment,
    }
    run_spec_digest = hashlib.sha256(json.dumps(
        run_spec, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")).hexdigest()
    artifact = {
        "artifact_schema_version": 1,
        "written_utc": datetime.now(timezone.utc).isoformat(),
        "run_spec_sha256": run_spec_digest,
        "provenance": run_spec,
        "result": dict(result),
    }
    serialized = json.dumps(artifact, sort_keys=True, indent=2, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(serialized)
        handle.flush()
        os.fsync(handle.fileno())
    reopened = json.loads(path.read_text(encoding="utf-8"))
    if reopened != artifact:
        raise RuntimeError("saved MU result did not match the in-memory result artifact")
    artifact_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "path": str(path),
        "artifact_sha256": artifact_hash,
        "run_spec_sha256": run_spec_digest,
        "size_bytes": path.stat().st_size,
    }


def _select_analysis_records(records, *, gap_sensitivity_only: bool):
    records = tuple(records)
    if not gap_sensitivity_only:
        return records
    return tuple(
        record for record in records
        if 78.0 < float(record.days[2] - record.days[1]) < 102.0
    )


def _validate_boundary_covariance(value: Mapping, *, name: str) -> dict[str, float]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must provide sigma_b_mm, tau_b_mm, and ell_mm")
    try:
        result = {key: float(value[key]) for key in ("sigma_b_mm", "tau_b_mm", "ell_mm")}
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must provide finite sigma_b_mm, tau_b_mm, and ell_mm") from exc
    if (not np.isfinite(tuple(result.values())).all()
            or result["sigma_b_mm"] < 0 or result["tau_b_mm"] < 0 or result["ell_mm"] <= 0):
        raise ValueError(f"{name} has invalid covariance values")
    if result["sigma_b_mm"] == 0 and result["tau_b_mm"] == 0:
        raise ValueError(f"{name} must include nonzero initial boundary variance")
    return result


def _volume_only_nuisance(records) -> dict:
    fitted = fit_mu_training_volume_nuisance(records)
    return {**fitted, "model_ready": True}


def _patient_bootstrap_mean_interval(values, *, seed: int, n_bootstrap: int = 10_000) -> dict:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        return {"n_patients": int(len(values)), "mean": None, "interval_95": None,
                "status": "fewer than two finite patient differences"}
    rng = np.random.default_rng(int(seed))
    indices = rng.integers(0, len(values), size=(int(n_bootstrap), len(values)))
    means = values[indices].mean(axis=1)
    lower, upper = np.quantile(means, [0.025, 0.975])
    return {"n_patients": int(len(values)), "mean": float(values.mean()),
            "interval_95": [float(lower), float(upper)], "n_bootstrap": int(n_bootstrap),
            "seed": int(seed), "bootstrap_unit": "patient"}


def _summarize_posterior_parameters(posterior: Mapping) -> dict:
    """Return compact history-only joint-parameter diagnostics, without posterior arrays."""
    grid = posterior.get("grid")
    points = np.asarray(grid.get("points") if isinstance(grid, Mapping) else None, dtype=float)
    weights = np.asarray(posterior.get("posterior_probabilities"), dtype=float)
    if (points.ndim != 2 or points.shape[1] != 2 or weights.shape != (len(points),)
            or len(points) < 2 or not np.isfinite(points).all() or not np.isfinite(weights).all()
            or np.any(weights < 0) or weights.sum() <= 0):
        raise ValueError("posterior diagnostics require finite joint log-parameter points and probabilities")
    weights = weights / weights.sum()
    mean = np.average(points, axis=0, weights=weights)
    centered = points - mean
    covariance = (centered * weights[:, None]).T @ centered
    standard_deviations = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    if np.all(standard_deviations > 0):
        correlation = float(covariance[0, 1] / np.prod(standard_deviations))
    else:
        correlation = None

    def weighted_quantiles(values):
        order = np.argsort(values, kind="mergesort")
        cumulative = np.cumsum(weights[order])
        return np.interp((0.025, 0.5, 0.975), cumulative, values[order])

    rho_quantiles = np.exp(weighted_quantiles(points[:, 0]))
    diffusion_quantiles = np.exp(weighted_quantiles(points[:, 1]))
    positive = weights > 0
    entropy = -float(np.sum(weights[positive] * np.log(weights[positive])))
    bound_mass = posterior.get("posterior_mass_at_parameter_bounds", {})
    if not isinstance(bound_mass, Mapping) or any(
        not np.isfinite(float(value)) or not 0 <= float(value) <= 1 for value in bound_mass.values()
    ):
        raise ValueError("posterior parameter-bound masses must be finite probabilities")
    return {
        "status": "descriptive posterior diagnostics; parameter identifiability is separate from forecast accuracy",
        "parameter_order": ["log_rho", "log_D"],
        "posterior_mean_log_parameters": [float(value) for value in mean],
        "posterior_covariance_log_parameters": covariance.astype(float).tolist(),
        "posterior_correlation_log_rho_log_D": correlation,
        "posterior_marginal_quantiles": {
            "rho_per_day_2_5_50_97_5": [float(value) for value in rho_quantiles],
            "D_mm2_per_day_2_5_50_97_5": [float(value) for value in diffusion_quantiles],
        },
        "posterior_effective_grid_cells": float(np.exp(entropy)),
        "posterior_mass_at_parameter_bounds": {
            str(key): float(value) for key, value in bound_mass.items()
        },
        "coarse_grid_points": int(posterior["coarse_grid_points"]),
        "refined_grid_points": int(posterior["refined_grid_points"]),
        "n_likelihood_evaluations": int(posterior["n_likelihood_evaluations"]),
        "posterior_draw_count": int(len(posterior.get("parameter_draws_log", ()))),
    }


def _hybrid_forecast(history: PatientHistory, nuisance: Mapping, domain_by_patient: Mapping,
                     boundary_covariance: Mapping, *, seed: int, settings: Mapping,
                     calibration_only: bool = False) -> dict:
    if not nuisance.get("model_ready", False):
        raise ForecastFailure(str(nuisance.get("unavailable_reason", "training nuisance fit unavailable")))
    domain = np.asarray(domain_by_patient[history.patient_id], dtype=bool)
    try:
        result = hybrid_posterior_predictive_forecast(
            masks=history.masks,
            days=history.days,
            target_day=history.forecast_day,
            prior={"distribution": "log_uniform"},
            h=2.0,
            boundary_covariance=dict(boundary_covariance),
            population_log_volume_slope_per_day=nuisance["population_log_volume_slope_per_day"],
            shrinkage_time_constant_days=nuisance["shrinkage_time_constant_days"],
            volume_variance_fit=nuisance["volume_variance_fit"],
            domain=domain,
            rho_bounds=tuple(settings["rho_bounds"]),
            D_bounds=tuple(settings["D_bounds"]),
            n_coarse=int(settings["n_coarse"]),
            mass_fraction=float(settings["mass_fraction"]),
            n_refine=int(settings["n_refine"]),
            n_draws=int(settings["n_draws"]),
            n_anchors=int(settings["n_anchors"]),
            seed=_stable_seed(seed, history.patient_id),
            calibration_only=calibration_only,
        )
    except ValueError as exc:
        if str(exc) in _EXPECTED_MODEL_FAILURES:
            raise ForecastFailure(str(exc)) from exc
        raise
    return result


def run_mu_first_two_to_third_development_cv(
    cohort, *, initial_boundary_covariance: Mapping,
    external_eligibility: ExternalEligibilityEvidence | None = None,
    gap_sensitivity_only: bool = False,
    seed: int = 20261028,
    boundary_calibration_seed: int = 20261029,
    calibration_minimum_patients: int = 8,
    calibration_folds: int = 5,
    rho_bounds=(1e-5, 1e-1), D_bounds=(1e-5, 10.0),
    n_coarse: int = 21, mass_fraction: float = 0.9, n_refine: int = 11,
    n_draws: int = 200, n_anchors: int = 256,
) -> dict:
    """Run one paired day-90 forecast per selected patient, with patient-level outer folds.

    ``cohort`` must come from ``load_mu_first_two_to_third``. Do not create it
    from masks directly: the loader enforces the source-backed evidence gate.
    The caller must prespecify ``initial_boundary_covariance`` before scoring;
    it seeds inner cross-fitted residual calibration and is not learned from a
    synthetic oracle. Patients with empty histories or unusable forecasts stay
    in the evaluation and receive the prespecified persistence fallback.
    Optional precision uses only a validated metadata-only external eligibility
    report; the runner does not accept a caller-entered external sample size.
    """
    required = ("records", "fold_by_patient", "origin_domain_by_patient",
                "target_valid_support_by_patient", "audit")
    if (not isinstance(cohort, MUDevelopmentCohort) or not is_loader_issued_cohort(cohort)
            or any(not hasattr(cohort, name) for name in required)):
        raise ValueError("cohort must be the exact in-memory object returned by the gated MU development loader")
    if cohort.audit.get("access", {}).get("source_backed_registration_gate_passed") is not True:
        raise PermissionError("scoring refused: source-backed registration evidence gate is not recorded as passed")
    if cohort.audit.get("access", {}).get("source_backed_target_support_gate_passed") is not True:
        raise PermissionError("scoring refused: source-backed target-support evidence gate is not recorded as passed")
    all_records = tuple(cohort.records)
    all_patient_ids = [str(record.patient_id) for record in all_records]
    if len(all_patient_ids) != 18 or len(set(all_patient_ids)) != len(all_patient_ids):
        raise ValueError("runner requires exactly the reviewed 18-patient MU cohort with unique IDs")
    if (set(cohort.fold_by_patient) != set(all_patient_ids)
            or set(cohort.origin_domain_by_patient) != set(all_patient_ids)
            or set(cohort.target_valid_support_by_patient) != set(all_patient_ids)):
        raise ValueError("folds and origin-only domains must cover exactly the loaded cohort")
    if not isinstance(gap_sensitivity_only, bool):
        raise ValueError("gap_sensitivity_only must be a boolean")
    if external_eligibility is not None and not is_validated_external_eligibility(external_eligibility):
        raise PermissionError("external precision requires a validated metadata-only eligibility report")
    if gap_sensitivity_only and external_eligibility is not None:
        raise ValueError("external precision is computed only for the full primary development cohort")
    if gap_sensitivity_only:
        records = _select_analysis_records(all_records, gap_sensitivity_only=True)
        analysis_subset = "central-gap sensitivity; excludes observed gaps 76-78 and 102-104 days inclusive"
    else:
        records = all_records
        analysis_subset = "all frozen 18 patients with observed gaps 76-104 days inclusive"
    patient_ids = [str(record.patient_id) for record in records]
    if len(patient_ids) < 5 or len(set(patient_ids)) != len(patient_ids):
        raise ValueError("analysis subset must contain at least five unique patients")
    analysis_folds = {patient_id: int(cohort.fold_by_patient[patient_id]) for patient_id in patient_ids}
    if set(analysis_folds.values()) != set(range(5)):
        raise ValueError("analysis subset must preserve at least one patient in every existing fold 0-4")
    analysis_domains = {patient_id: cohort.origin_domain_by_patient[patient_id] for patient_id in patient_ids}
    analysis_support = {patient_id: cohort.target_valid_support_by_patient[patient_id] for patient_id in patient_ids}
    if (isinstance(calibration_minimum_patients, bool)
            or int(calibration_minimum_patients) != calibration_minimum_patients
            or int(calibration_minimum_patients) < 5):
        raise ValueError("calibration_minimum_patients must be at least five")
    if (isinstance(calibration_folds, bool) or int(calibration_folds) != calibration_folds
            or int(calibration_folds) < 2):
        raise ValueError("calibration_folds must be at least two")
    initial_covariance = _validate_boundary_covariance(
        initial_boundary_covariance, name="initial_boundary_covariance"
    )
    settings = {
        "rho_bounds": tuple(rho_bounds), "D_bounds": tuple(D_bounds),
        "n_coarse": int(n_coarse), "mass_fraction": float(mass_fraction),
        "n_refine": int(n_refine), "n_draws": int(n_draws), "n_anchors": int(n_anchors),
    }
    if (settings["n_coarse"] < 2 or settings["n_refine"] < 2 or settings["n_draws"] < 2
            or settings["n_anchors"] < 3 or not 0 < settings["mass_fraction"] < 1):
        raise ValueError("forecast grid, draw, anchor, or posterior-mass settings are invalid")

    fold_fit_audits = []

    def fit_outer_nuisance(training_records):
        training_records = tuple(training_records)
        try:
            volume_nuisance = _volume_only_nuisance(training_records)
        except ValueError as exc:
            expected = (
                "need at least three positive-volume training patients",
                "training records must provide at least two effective forecast horizons",
            )
            if not any(message in str(exc) for message in expected):
                raise
            fold_fit_audits.append({
                "training_patient_count": len(training_records),
                "model_ready": False,
                "unavailable_reason": str(exc),
                "boundary_calibration_patient_count": 0,
            })
            return {"model_ready": False, "unavailable_reason": str(exc)}

        # Boundary residuals need a genuine surface. This rule affects only
        # calibration inside the outer-training fold; outer test rows remain.
        calibration_records = tuple(
            record for record in training_records
            if all(np.any(mask) and not np.all(mask) for mask in record.masks)
        )
        if len(calibration_records) < int(calibration_minimum_patients):
            reason = "too few training patients with three usable nonempty, nonfull masks for boundary calibration"
            fold_fit_audits.append({
                "training_patient_count": len(training_records),
                "model_ready": False,
                "unavailable_reason": reason,
                "boundary_calibration_patient_count": len(calibration_records),
            })
            return {"model_ready": False, "unavailable_reason": reason}

        inner_folds = min(int(calibration_folds), len(calibration_records))

        def inner_forecast(history, inner_nuisance):
            return _hybrid_forecast(
                history, inner_nuisance, analysis_domains,
                initial_covariance, seed=int(boundary_calibration_seed), settings=settings,
                calibration_only=True,
            )

        def inner_fallback(history, _inner_nuisance):
            return {"pre_boundary_probability": np.asarray(history.masks[-1], dtype=np.float32),
                    "_calibration_fallback": True}

        def inner_nuisance_fit(inner_training_records):
            return _volume_only_nuisance(inner_training_records)

        def inner_forecast_with_patient_capture(history, inner_nuisance):
            prediction = inner_forecast(history, inner_nuisance)
            return prediction

        # The calibration helper owns residual construction. Failed predictions
        # use persistence but are marked and excluded from covariance fitting.
        try:
            calibration = run_cross_fitted_boundary_calibration(
                calibration_records,
                fit_nuisance=inner_nuisance_fit,
                forecast=inner_forecast_with_patient_capture,
                n_folds=inner_folds,
                seed=int(boundary_calibration_seed),
                minimum_history=2,
                spacing_mm=2.0,
                n_anchors=int(settings["n_anchors"]),
                ell_candidates_mm=(5.0, 10.0, 20.0),
                forecast_horizon_days=90.0,
                fallback_forecast=inner_fallback,
                allow_empty_history_masks=True,
            )
        except ValueError as exc:
            expected = (
                "need at least two patient-held-out calibration cases",
                "need calibration origins from at least two distinct patients",
                "need at least three positive-volume training patients",
                "training records must provide at least two effective forecast horizons",
            )
            if not any(message in str(exc) for message in expected):
                raise
            reason = f"cross-fitted boundary calibration unavailable: {exc}"
            fold_fit_audits.append({
                "training_patient_count": len(training_records),
                "model_ready": False,
                "unavailable_reason": reason,
                "boundary_calibration_patient_count": len(calibration_records),
            })
            return {"model_ready": False, "unavailable_reason": reason}
        inner_fallback_count = sum(row["fallback_count"] for row in calibration["fold_audit"])
        if inner_fallback_count:
            reason = "inner cross-fitted calibration used persistence fallback; boundary estimate withheld"
            fold_fit_audits.append({
                "training_patient_count": len(training_records),
                "model_ready": False,
                "unavailable_reason": reason,
                "boundary_calibration_patient_count": len(calibration_records),
                "valid_boundary_residual_count": calibration["n_residual_cases"],
                "calibration_fallback_count": inner_fallback_count,
            })
            return {"model_ready": False, "unavailable_reason": reason}
        covariance = calibration["covariance_fit"]
        if calibration["n_residual_cases"] < 2:
            reason = "cross-fitted boundary calibration produced fewer than two valid patient residuals"
            fold_fit_audits.append({
                "training_patient_count": len(training_records),
                "model_ready": False,
                "unavailable_reason": reason,
                "boundary_calibration_patient_count": len(calibration_records),
                "valid_boundary_residual_count": calibration["n_residual_cases"],
            })
            return {"model_ready": False, "unavailable_reason": reason}
        fitted = {
            **volume_nuisance,
            "model_ready": True,
            "boundary_covariance": {
                "sigma_b_mm": covariance["sigma_b_mm"],
                "tau_b_mm": covariance["tau_b_mm"],
                "ell_mm": covariance["ell_mm"],
            },
        }
        fold_fit_audits.append({
            "training_patient_count": len(training_records),
            "model_ready": True,
            "volume_nuisance_audit": volume_nuisance["audit_counts"],
            "boundary_calibration_patient_count": len(calibration_records),
            "valid_boundary_residual_count": calibration["n_residual_cases"],
            "boundary_covariance": fitted["boundary_covariance"],
            "boundary_covariance_candidates": covariance["candidate_fits"],
            "boundary_length_scale_on_search_boundary": covariance["selection_boundary"],
            "boundary_calibration_scope": calibration["scope"],
            "calibration_fallback_count": inner_fallback_count,
        })
        return fitted

    parameter_diagnostics_by_patient = {}

    def forecast(history, nuisance):
        if not nuisance.get("model_ready", False):
            raise ForecastFailure(str(nuisance.get("unavailable_reason", "outer-fold model unavailable")))
        if any(not np.asarray(mask, dtype=bool).any() for mask in history.masks):
            raise ForecastFailure("empty history mask requires the prespecified persistence fallback")
        predicted = _hybrid_forecast(
            history, nuisance, analysis_domains,
            nuisance["boundary_covariance"], seed=int(seed), settings=settings,
        )
        parameter_diagnostics_by_patient[history.patient_id] = _summarize_posterior_parameters(
            predicted["posterior"]
        )
        return {"forecast": predicted, "history_mask": history.masks[-1]}

    def persistence_fallback(history, _nuisance):
        last = np.asarray(history.masks[-1], dtype=bool)
        return {"fallback_masks": {"hybrid": last, "no_pde": last, "persistence": last,
                                   "geometric_train_rate": last},
                "history_mask": last}

    cv = run_first_third_90day_patient_cv(
        records, analysis_folds, fit_outer_nuisance, forecast,
        make_first_third_score_callback(
            spacing_mm=2.0,
            valid_support_by_patient=analysis_support,
        ),
        fallback_forecast=persistence_fallback,
        n_folds=5,
    )
    score_by_id = {row["patient_id"]: row["scores"] for row in cv["rows"]}
    if set(score_by_id) != set(patient_ids) or len(score_by_id) != len(patient_ids):
        raise RuntimeError("CV output did not preserve exactly one row per selected patient")
    differences = [(patient_id, score_by_id[patient_id]["paired_delta_dice_hybrid_minus_no_pde"])
                   for patient_id in patient_ids]
    gap_by_id = {record.patient_id: float(record.days[2] - record.days[1]) for record in records}
    for row in cv["rows"]:
        row["observed_target_gap_days"] = gap_by_id[row["patient_id"]]
        row["parameter_posterior_diagnostics"] = parameter_diagnostics_by_patient.get(
            row["patient_id"], {"status": "unavailable; forecast used persistence fallback"},
        )
    primary_delta_values = [value for _, value in differences]
    precision = None
    if external_eligibility is not None:
        precision = precision_gate_from_patient_differences(
            differences,
            expected_patient_ids=patient_ids,
            external_eligibility=external_eligibility,
        )

    return {
        "schema_version": 1,
        "study": "exploratory MU first-two-to-third PDE development",
        "protocol": {
            "analysis_subset": analysis_subset,
            "forecast_horizon_days": 90.0,
            "observed_followup_window_days": [76.0, 104.0],
            "primary_difference": "Dice(hybrid) - Dice(matched-volume no-PDE)",
            "comparator": "signed-distance volume-ranker with matched volume draws and same boundary covariance",
            "fallback": "persistence for hybrid, no-PDE, and persistence arms; paired difference zero",
            "log_parameter_prior": "fixed log-uniform",
            "initial_boundary_covariance": initial_covariance,
            "initial_boundary_covariance_role": "prespecified seed for nested training-only boundary calibration; not a fitted truth",
            "spatial_resolution_mm": 2.0,
            "ensemble_inclusion_threshold": 0.5,
            "empty_mask_dice_rule": "both_empty=1; one_empty=0",
            "forecast_settings": settings,
        },
        "random_seeds": {
            "forecast_base": int(seed),
            "boundary_calibration": int(boundary_calibration_seed),
            "mu_paired_difference_bootstrap": 20261030,
            "geometric_train_rate_bootstrap": 20261031,
            "external_precision_bootstrap": 20261028,
        },
        "cohort_audit": cohort.audit,
        "outer_fold_fit_audits": fold_fit_audits,
        "fold_audit": cv["fold_audit"],
        "rows": cv["rows"],
        "n_patients": len(cv["rows"]),
        "n_successful_hybrid_forecasts": sum(row["forecast_status"] == "ok" for row in cv["rows"]),
        "n_persistence_fallbacks": sum(row["forecast_status"] == "fallback" for row in cv["rows"]),
        "n_geometric_train_rate_persistence_fallbacks": sum(
            row["scores"]["geometric_train_rate_status"] != "forecast" for row in cv["rows"]
        ),
        "hybrid_forecast_success_rate": float(
            sum(row["forecast_status"] == "ok" for row in cv["rows"]) / len(cv["rows"])
        ),
        "n_empty_targets": sum(bool(row["scores"]["target_empty"]) for row in cv["rows"]),
        "primary_paired_difference_bootstrap": _patient_bootstrap_mean_interval(
            primary_delta_values, seed=20261030,
        ),
        "geometric_train_rate_paired_difference_bootstrap": _patient_bootstrap_mean_interval(
            [row["scores"]["paired_delta_dice_hybrid_minus_geometric_train_rate"] for row in cv["rows"]],
            seed=20261031,
        ),
        "boundary_gap_sensitivity_status": (
            "this run used a central-gap-only cohort and refit all outer/inner nuisances on that subset"
            if gap_sensitivity_only else
            "not run here; use the protocol-complete wrapper for the central-gap nested refit"
        ),
        "prior_support_sensitivity_status": (
            "not run by this primary runner; freeze and run one-decade narrower and wider bounds before scoring"
        ),
        "precision_gate": precision,
        "precision_gate_status": (
            "not computed: validated metadata-only external eligibility report unavailable"
            if external_eligibility is None else precision["interpretation"]
        ),
        "interpretation_limit": (
            "Previously explored MU data provide exploratory development evidence only. This CV does not establish "
            "independent confirmation, external performance, clinical validity, or novelty. A boundary covariance "
            "search-boundary result is inconclusive; report every fallback and every patient-level paired difference."
        ),
    }


def run_mu_first_two_to_third_with_prior_sensitivity(
    cohort, *,
    external_eligibility: ExternalEligibilityEvidence | None = None,
    **runner_options,
) -> dict:
    """Run three prior-support evaluations plus a central-gap-only nested refit.

    This is the protocol-complete entry point. It runs four full nested patient-CV analyses. The three prior-support
    runs use the same full cohort; the gap sensitivity excludes edge-gap patients from every outer and inner fit.
    All runs preserve the existing fold assignments and random seeds. Only the primary bounds receive the optional
    validated external-eligibility evidence for precision calculation. This can be substantially slower than one
    primary run.
    """
    if ("rho_bounds" in runner_options or "D_bounds" in runner_options
            or "gap_sensitivity_only" in runner_options
            or "initial_boundary_covariance" in runner_options):
        raise ValueError(
            "the protocol-complete wrapper owns prior bounds, boundary seed, and gap-sensitivity selection"
        )
    scenarios = {
        "primary": ((1e-5, 1e-1), (1e-5, 10.0)),
        "one_decade_narrower": ((1e-4, 1e-2), (1e-4, 1.0)),
        "one_decade_wider": ((1e-6, 1.0), (1e-6, 100.0)),
    }
    runs = {}
    for name, (rho_bounds, D_bounds) in scenarios.items():
        runs[name] = run_mu_first_two_to_third_development_cv(
            cohort,
            initial_boundary_covariance=LOCKED_MU_INITIAL_BOUNDARY_COVARIANCE,
            external_eligibility=(external_eligibility if name == "primary" else None),
            rho_bounds=rho_bounds,
            D_bounds=D_bounds,
            **runner_options,
        )

    central_gap_run = run_mu_first_two_to_third_development_cv(
        cohort,
        initial_boundary_covariance=LOCKED_MU_INITIAL_BOUNDARY_COVARIANCE,
        external_eligibility=None,
        gap_sensitivity_only=True,
        rho_bounds=scenarios["primary"][0],
        D_bounds=scenarios["primary"][1],
        **runner_options,
    )

    patient_sets = [{row["patient_id"] for row in run["rows"]} for run in runs.values()]
    if any(patient_set != patient_sets[0] for patient_set in patient_sets[1:]):
        raise RuntimeError("prior-support sensitivity runs did not preserve the same patient set")
    primary_delta = {row["patient_id"]: row["scores"]["paired_delta_dice_hybrid_minus_no_pde"]
                     for row in runs["primary"]["rows"]}
    sensitivity = {}
    for name, run in runs.items():
        deltas = {row["patient_id"]: row["scores"]["paired_delta_dice_hybrid_minus_no_pde"]
                  for row in run["rows"]}
        sensitivity[name] = {
            "rho_bounds": list(scenarios[name][0]),
            "D_bounds": list(scenarios[name][1]),
            "mean_paired_delta_dice": float(np.mean(list(deltas.values()))),
            "n_successful_hybrid_forecasts": run["n_successful_hybrid_forecasts"],
            "n_persistence_fallbacks": run["n_persistence_fallbacks"],
            "mean_delta_change_from_primary": float(np.mean([
                deltas[patient_id] - primary_delta[patient_id] for patient_id in sorted(primary_delta)
            ])),
            "direction_matches_primary": bool(
                np.sign(np.mean(list(deltas.values()))) == np.sign(np.mean(list(primary_delta.values())))
            ),
        }
    return {
        "schema_version": 1,
        "study": "exploratory MU first-two-to-third PDE development with prior-support sensitivity",
        "patient_set_sha256": runs["primary"]["cohort_audit"]["source_hashes"][
            "selected_patient_ids_sha256_sorted_nul_joined"
        ],
        "runs": runs,
        "prior_support_sensitivity": sensitivity,
        "prior_support_sensitivity_status": "completed across identical patient folds and random seeds",
        "boundary_gap_sensitivity_run": central_gap_run,
        "boundary_gap_sensitivity_status": (
            "completed as a central-gap-only nested refit; edge-gap patients were excluded from outer and inner fitting"
        ),
        "interpretation_limit": (
            "All runs use previously explored MU development data. They do not establish independent confirmation, "
            "external performance, clinical validity, or novelty."
        ),
    }


def run_and_write_mu_first_two_to_third_with_prior_sensitivity(
    cohort, *, output_path: str | Path,
    external_eligibility: ExternalEligibilityEvidence | None = None,
    **runner_options,
) -> dict:
    """Run the locked four-analysis MU development workflow and save one immutable result artifact."""
    result = run_mu_first_two_to_third_with_prior_sensitivity(
        cohort,
        external_eligibility=external_eligibility,
        **runner_options,
    )
    artifact = write_mu_development_result(result, output_path)
    return {"result": result, "artifact_write": artifact}
