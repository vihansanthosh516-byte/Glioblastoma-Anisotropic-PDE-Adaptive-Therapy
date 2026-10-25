#!/usr/bin/env python3
"""Exploratory patient-grouped CV for the Track B growth-zone hybrid.

This runner is intentionally separate from the existing whole-mask pipeline.
It requires an explicit acknowledgement flag before loading/scoring patient
images and writes no voxel-level predictions.  All MU results are marked
exploratory because this cohort has already been examined in earlier Track B
work; a new fold assignment cannot make it an independent confirmation set.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold


ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = ROOT / "data/tcia/MU-Glioma-Post"
PATIENT_MANIFEST = ROOT / "data/manifests/patient_manifest.csv"
PAIR_MANIFEST = ROOT / "data/manifests/forecast_pairs_mu.csv"


def _load_local_module(alias: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(alias, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load local module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[alias] = module
    spec.loader.exec_module(module)
    return module


FEATURES = _load_local_module("growth_zone_features_impl", ROOT / "src/160_growth_zone_features.py")
MODELS = _load_local_module("growth_zone_model_impl", ROOT / "src/161_growth_zone_model.py")


@dataclass(frozen=True)
class RawSample:
    """Forecast features and target, before fold-specific gate scaling."""

    sample_id: str
    patient_id: str
    features: np.ndarray
    log_pde: np.ndarray
    target: np.ndarray
    log_instability: np.ndarray
    prior_instability: np.ndarray
    prior_addition_count: int
    prior_candidate_voxels: int
    prior_interval_days: float
    forecast_interval_days: float
    candidate_count: int
    all_addition_count: int
    outside_roi_addition_count: int

    def materialize(self, scales: tuple[float, float], strength: float = 1.0) -> Any:
        log_scale, prior_scale = scales
        instability = 0.5 * np.clip(self.log_instability / log_scale, 0.0, 1.0)
        instability += 0.5 * np.clip(self.prior_instability / prior_scale, 0.0, 1.0)
        gate = np.exp(-strength * instability).astype(np.float32, copy=False)
        return MODELS.VoxelSample(
            patient_id=self.patient_id,
            features=self.features,
            log_pde=self.log_pde,
            gate=gate,
            target=self.target.astype(np.uint8, copy=False),
            prior_addition_count=self.prior_addition_count,
            prior_candidate_voxels=self.prior_candidate_voxels,
            prior_interval_days=self.prior_interval_days,
            forecast_interval_days=self.forecast_interval_days,
        )


def _truth(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _manifest_triplets(patient_manifest: Path, pair_manifest: Path) -> list[tuple[str, dict[str, str], dict[str, str]]]:
    patients = {
        row["patient_id"]: row
        for row in _read_csv(patient_manifest)
        if _truth(row.get("eligible", "")) and _truth(row.get("is_gbm", ""))
    }
    by_patient: dict[str, list[dict[str, str]]] = {}
    for row in _read_csv(pair_manifest):
        patient_id = row.get("patient_id", "")
        if patient_id not in patients or not _truth(row.get("in_primary", "")):
            continue
        if not _truth(row.get("input_core_nonempty", "")):
            continue
        by_patient.setdefault(patient_id, []).append(row)

    triplets: list[tuple[str, dict[str, str], dict[str, str]]] = []
    for patient_id, pairs in by_patient.items():
        pairs.sort(key=lambda row: (float(row["day_in"]), int(row["tp_in"]), int(row["tp_out"])))
        for first, second in zip(pairs, pairs[1:]):
            if int(first["tp_out"]) != int(second["tp_in"]):
                continue
            if abs(float(first["day_out"]) - float(second["day_in"])) > 1e-3:
                continue
            if not _truth(second.get("input_core_nonempty", "")):
                continue
            triplets.append((patient_id, first, second))
    if not triplets:
        raise ValueError("no chained, eligible GBM scan triplets found in the manifests")
    return triplets


def _single_file(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise ValueError(f"expected one file matching {pattern!r} in {directory}, found {len(matches)}")
    return matches[0]


def _load_mask(data_root: Path, patient_id: str, timepoint: int) -> tuple[np.ndarray, np.ndarray]:
    import nibabel as nib

    directory = data_root / patient_id / f"Timepoint_{timepoint}"
    path = _single_file(directory, "*_tumorMask.nii.gz")
    image = nib.load(str(path))
    labels = np.rint(np.asanyarray(image.dataobj)).astype(np.int16, copy=False)
    return labels, np.asarray(image.affine, dtype=np.float64)


def _load_visit(data_root: Path, patient_id: str, timepoint: int) -> Any:
    import nibabel as nib

    directory = data_root / patient_id / f"Timepoint_{timepoint}"
    mask_path = _single_file(directory, "*_tumorMask.nii.gz")
    mask_image = nib.load(str(mask_path))
    labels = np.rint(np.asanyarray(mask_image.dataobj)).astype(np.int16, copy=False)
    affine = np.asarray(mask_image.affine, dtype=np.float64)
    images: dict[str, np.ndarray] = {}
    for modality in FEATURES.MODALITIES:
        path = _single_file(directory, f"*_brain_{modality}.nii.gz")
        image = nib.load(str(path))
        if image.shape != labels.shape or not np.allclose(image.affine, affine, rtol=0.0, atol=1e-3):
            raise ValueError(f"{patient_id} timepoint {timepoint} {modality} is not aligned to its mask")
        images[modality] = np.asarray(image.dataobj, dtype=np.float32)
    return FEATURES.VisitVolume(images=images, labels=labels, affine=affine)


def _load_raw_samples(
    triplets: Sequence[tuple[str, dict[str, str], dict[str, str]]],
    data_root: Path,
    diffusion: Any,
) -> list[RawSample]:
    samples: list[RawSample] = []
    for patient_id, first, second in triplets:
        tp1, tp2, tp3 = int(first["tp_in"]), int(first["tp_out"]), int(second["tp_out"])
        previous = _load_visit(data_root, patient_id, tp1)
        current = _load_visit(data_root, patient_id, tp2)
        future_labels, future_affine = _load_mask(data_root, patient_id, tp3)
        if future_labels.shape != current.labels.shape or not np.allclose(
            future_affine, current.affine, rtol=0.0, atol=1e-3
        ):
            raise ValueError(f"{patient_id} timepoint {tp3} is not aligned to the forecast origin")

        visit_features = FEATURES.build_visit_features(previous, current, diffusion)
        log_draws, prior_draws = FEATURES.fixed_perturbation_draws(
            previous, current, visit_features, diffusion
        )
        log_instability, prior_instability = FEATURES.perturbation_instability_maps(
            visit_features.log_pde_field,
            log_draws,
            visit_features.recent_additions,
            prior_draws,
            visit_features.candidate,
        )
        x, _ = FEATURES.feature_matrix(visit_features)
        target_full = FEATURES.make_transition_target(
            current.labels, future_labels, visit_features.candidate
        )
        target = target_full[visit_features.candidate]
        all_additions = (
            FEATURES.annotated_core_mask(future_labels)
            & ~FEATURES.annotated_core_mask(current.labels)
        )
        target_count = int(np.count_nonzero(target))
        all_count = int(np.count_nonzero(all_additions))
        sample_id = f"{patient_id}:{tp1}-{tp2}-{tp3}"
        samples.append(
            RawSample(
                sample_id=sample_id,
                patient_id=patient_id,
                features=x,
                log_pde=visit_features.log_pde_field[visit_features.candidate].astype(np.float32),
                target=target.astype(np.uint8),
                log_instability=log_instability[visit_features.candidate].astype(np.float32),
                prior_instability=prior_instability[visit_features.candidate].astype(np.float32),
                prior_addition_count=visit_features.prior_addition_count,
                prior_candidate_voxels=visit_features.prior_candidate_voxels,
                prior_interval_days=float(first["dt_days"]),
                forecast_interval_days=float(second["dt_days"]),
                candidate_count=int(np.count_nonzero(visit_features.candidate)),
                all_addition_count=all_count,
                outside_roi_addition_count=max(0, all_count - target_count),
            )
        )
    return samples


def _fold_scales(samples: Sequence[RawSample]) -> tuple[float, float]:
    by_patient: dict[str, list[float]] = {}
    by_patient_prior: dict[str, list[float]] = {}
    for sample in samples:
        by_patient.setdefault(sample.patient_id, []).append(
            float(np.percentile(sample.log_instability, 90.0)) if sample.log_instability.size else 0.0
        )
        by_patient_prior.setdefault(sample.patient_id, []).append(
            float(np.percentile(sample.prior_instability, 90.0)) if sample.prior_instability.size else 0.0
        )
    log_scale = float(np.median([np.median(values) for values in by_patient.values()]))
    prior_scale = float(np.median([np.median(values) for values in by_patient_prior.values()]))
    return max(log_scale, 1e-6), max(prior_scale, 1e-6)


def _mean_patient_metric(samples: Sequence[Any], values: Sequence[float]) -> float:
    grouped: dict[str, list[float]] = {}
    for sample, value in zip(samples, values):
        if np.isfinite(value):
            patient_id = sample["patient_id"] if isinstance(sample, dict) else sample.patient_id
            grouped.setdefault(patient_id, []).append(float(value))
    if not grouped:
        return float("nan")
    return float(np.mean([np.mean(patient_values) for patient_values in grouped.values()]))


def _sample_ap(sample: Any, score: np.ndarray, count: float) -> tuple[np.ndarray, dict[str, float]]:
    probability = MODELS.calibrate_probability_mass(score, count)
    target = np.asarray(sample.target, dtype=np.uint8)
    metrics = {
        "average_precision": float(average_precision_score(target, probability))
        if np.unique(target).size == 2 else float("nan"),
        "auc": float(roc_auc_score(target, probability)) if np.unique(target).size == 2 else float("nan"),
        "brier": float(brier_score_loss(target, probability)),
        "predicted_mass": float(probability.sum()),
    }
    return probability, metrics


def _choose_ridge(
    train_raw: Sequence[RawSample],
    ridge_grid: Sequence[float],
    inner_folds: int,
    seed: int,
) -> float:
    groups = np.asarray([sample.patient_id for sample in train_raw])
    n_groups = np.unique(groups).size
    n_folds = min(inner_folds, n_groups)
    if n_folds < 2:
        raise ValueError("too few training patients for inner grouped selection")
    splitter = GroupKFold(n_splits=n_folds)
    scores = {float(ridge): [] for ridge in ridge_grid}
    for inner_train_idx, inner_val_idx in splitter.split(np.zeros(len(train_raw)), groups=groups):
        inner_train_raw = [train_raw[index] for index in inner_train_idx]
        inner_val_raw = [train_raw[index] for index in inner_val_idx]
        scales = _fold_scales(inner_train_raw)
        inner_train = [sample.materialize(scales) for sample in inner_train_raw]
        inner_val = [sample.materialize(scales) for sample in inner_val_raw]
        prior = MODELS.fit_count_prior(inner_train)
        for ridge in ridge_grid:
            model = MODELS.fit_residual_model(
                inner_train,
                ridge=float(ridge),
                seed=seed,
                use_pde=True,
                use_gate=True,
            )
            ap = []
            for raw, sample in zip(inner_val_raw, inner_val):
                count = MODELS.predict_addition_count(
                    prior,
                    sample.prior_addition_count,
                    sample.prior_candidate_voxels,
                    sample.prior_interval_days,
                    sample.forecast_interval_days,
                    sample.features.shape[0],
                )
                score = model.predict_score(sample.features, sample.log_pde, sample.gate)
                _, metric = _sample_ap(sample, score, count)
                ap.append(metric["average_precision"])
            scores[float(ridge)].append(_mean_patient_metric(inner_val, ap))
    means = {
        ridge: float(np.nanmean(values)) if np.any(np.isfinite(values)) else float("-inf")
        for ridge, values in scores.items()
    }
    return max(means, key=means.get)


def _choose_tree_leaves(
    train_raw: Sequence[RawSample],
    inner_folds: int,
    seed: int,
    leaf_options: Sequence[int] = (7, 15),
) -> int:
    """Choose nonlinear-baseline complexity using training patients only."""
    groups = np.asarray([sample.patient_id for sample in train_raw])
    n_groups = np.unique(groups).size
    n_folds = min(inner_folds, n_groups)
    if n_folds < 2:
        raise ValueError("too few training patients for inner grouped tree selection")
    splitter = GroupKFold(n_splits=n_folds)
    scores = {int(leaves): [] for leaves in leaf_options}
    for inner_train_idx, inner_val_idx in splitter.split(np.zeros(len(train_raw)), groups=groups):
        inner_train_raw = [train_raw[index] for index in inner_train_idx]
        inner_val_raw = [train_raw[index] for index in inner_val_idx]
        scales = _fold_scales(inner_train_raw)
        inner_train = [sample.materialize(scales) for sample in inner_train_raw]
        inner_val = [sample.materialize(scales) for sample in inner_val_raw]
        count_prior = MODELS.fit_count_prior(inner_train)
        for leaves in leaf_options:
            model = MODELS.fit_tree_baseline(
                inner_train,
                max_leaf_nodes=int(leaves),
                max_per_class_per_interval=1000,
                seed=seed,
            )
            ap = []
            for raw, sample in zip(inner_val_raw, inner_val):
                count = MODELS.predict_addition_count(
                    count_prior,
                    sample.prior_addition_count,
                    sample.prior_candidate_voxels,
                    sample.prior_interval_days,
                    sample.forecast_interval_days,
                    sample.features.shape[0],
                )
                score = model.predict_score(sample.features)
                _, metric = _sample_ap(sample, score, count)
                ap.append(metric["average_precision"])
            scores[int(leaves)].append(_mean_patient_metric(inner_val, ap))
    means = {
        leaves: float(np.nanmean(values)) if np.any(np.isfinite(values)) else float("-inf")
        for leaves, values in scores.items()
    }
    return max(means, key=means.get)


def _evaluate_outer_fold(
    train_raw: Sequence[RawSample],
    test_raw: Sequence[RawSample],
    ridge_grid: Sequence[float],
    inner_folds: int,
    seed: int,
    fold: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    scales = _fold_scales(train_raw)
    train = [sample.materialize(scales) for sample in train_raw]
    test = [sample.materialize(scales) for sample in test_raw]
    ridge = _choose_ridge(train_raw, ridge_grid, inner_folds, seed + fold)
    tree_leaves = _choose_tree_leaves(train_raw, inner_folds, seed + fold)
    count_prior = MODELS.fit_count_prior(train)
    full_model = MODELS.fit_residual_model(train, ridge=ridge, seed=seed + fold)
    no_gate_model = MODELS.fit_residual_model(
        train, ridge=ridge, seed=seed + fold, use_pde=True, use_gate=False
    )
    matched_model = MODELS.fit_residual_model(
        train, ridge=ridge, seed=seed + fold, use_pde=False, use_gate=False
    )
    tree_model = MODELS.fit_tree_baseline(
        train, max_leaf_nodes=tree_leaves, max_per_class_per_interval=1000, seed=seed + fold
    )
    no_history_train = [replace(sample, features=np.column_stack((sample.features[:, :-1], np.zeros(len(sample.target))))) for sample in train]
    no_history_model = MODELS.fit_residual_model(no_history_train, ridge=ridge, seed=seed + fold)

    fold_rows: list[dict[str, Any]] = []
    for raw, sample in zip(test_raw, test):
        predicted_count = MODELS.predict_addition_count(
            count_prior,
            sample.prior_addition_count,
            sample.prior_candidate_voxels,
            sample.prior_interval_days,
            sample.forecast_interval_days,
            raw.candidate_count,
        )
        cohort_count = MODELS.predict_cohort_rate_count(
            count_prior, sample.forecast_interval_days, raw.candidate_count
        )
        hybrid_score = full_model.predict_score(sample.features, sample.log_pde, sample.gate)
        no_gate_score = no_gate_model.predict_score(sample.features, sample.log_pde, sample.gate)
        matched_score = matched_model.predict_score(sample.features, sample.log_pde, sample.gate)
        tree_score = tree_model.predict_score(sample.features)
        no_history_features = np.column_stack((sample.features[:, :-1], np.zeros(len(sample.target))))
        no_history_score = no_history_model.predict_score(no_history_features, sample.log_pde, sample.gate)
        distance_score = -sample.features[:, 0]
        spatial_scores = {
            "hybrid": hybrid_score,
            "hybrid_no_stability_gate": no_gate_score,
            "matched_logistic_no_pde": matched_score,
            "hist_gradient_boosting_no_pde": tree_score,
            "hybrid_no_prior_transition": no_history_score,
            "distance_only": distance_score,
            "pde_only": sample.log_pde,
        }
        arm_metrics: dict[str, Any] = {}
        cohort_count_brier = None
        for arm, score in spatial_scores.items():
            _, arm_metrics[arm] = _sample_ap(sample, score, predicted_count)
            if arm == "hybrid":
                _, cohort_metrics = _sample_ap(sample, score, cohort_count)
                cohort_count_brier = cohort_metrics["brier"]
        target_count = int(np.sum(sample.target))
        fold_rows.append(
            {
                "patient_id": raw.patient_id,
                "sample_id": raw.sample_id,
                "fold": fold,
                "n_candidates": raw.candidate_count,
                "target_additions_in_roi": target_count,
                "target_additions_total": raw.all_addition_count,
                "additions_outside_roi": raw.outside_roi_addition_count,
                "predicted_addition_count": predicted_count,
                "cohort_rate_count": cohort_count,
                "count_absolute_error": abs(predicted_count - target_count),
                "cohort_count_absolute_error": abs(cohort_count - target_count),
                "hybrid_cohort_count_brier": cohort_count_brier,
                "arms": arm_metrics,
            }
        )
    return fold_rows, {"fold": fold, "selected_ridge": ridge, "selected_tree_max_leaf_nodes": tree_leaves, "n_train_patients": len({s.patient_id for s in train_raw}), "n_test_patients": len({s.patient_id for s in test_raw}), "gate_scales": {"pde": scales[0], "history": scales[1]}}


def _paired_patient_bootstrap(
    rows: Sequence[dict[str, Any]],
    arm_a: str,
    arm_b: str,
    metric: str,
    seed: int,
    replicates: int = 2000,
) -> dict[str, float | int]:
    by_patient: dict[str, dict[str, list[float]]] = {}
    for row in rows:
        entry = by_patient.setdefault(row["patient_id"], {arm_a: [], arm_b: []})
        entry[arm_a].append(row["arms"][arm_a][metric])
        entry[arm_b].append(row["arms"][arm_b][metric])
    differences = []
    for entry in by_patient.values():
        a = np.asarray(entry[arm_a], dtype=float)
        b = np.asarray(entry[arm_b], dtype=float)
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if a.size and b.size:
            differences.append(float(np.mean(a) - np.mean(b)))
    if len(differences) < 2:
        return {"n_patients": len(differences), "mean_difference": float("nan"), "ci95_low": float("nan"), "ci95_high": float("nan")}
    values = np.asarray(differences, dtype=float)
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(replicates, values.size), replace=True).mean(axis=1)
    return {
        "n_patients": int(values.size),
        "mean_difference": float(values.mean()),
        "ci95_low": float(np.percentile(draws, 2.5)),
        "ci95_high": float(np.percentile(draws, 97.5)),
    }


def _json_safe(value: Any) -> Any:
    """Convert NumPy values and undefined metrics to strict JSON values."""
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def run(args: argparse.Namespace) -> dict[str, Any]:
    triplets = _manifest_triplets(args.patient_manifest, args.pair_manifest)
    diffusion = FEATURES.DiffusionParameters(
        edge_sensitivity=args.edge_sensitivity,
        screening_length_mm=args.screening_length_mm,
        label2_log_permeability=args.label2_log_permeability,
        candidate_radius_mm=args.candidate_radius_mm,
        prior_smoothing_mm=args.prior_smoothing_mm,
        block_cavity=not args.open_cavity_ablation,
    )
    raw_samples = _load_raw_samples(triplets, args.data_root, diffusion)
    if len({sample.patient_id for sample in raw_samples}) < args.outer_folds:
        raise ValueError("fewer unique patients than outer folds")
    groups = np.asarray([sample.patient_id for sample in raw_samples])
    splitter = GroupKFold(n_splits=args.outer_folds)
    ridge_grid = tuple(float(value) for value in args.ridge_grid)
    rows: list[dict[str, Any]] = []
    folds: list[dict[str, Any]] = []
    for fold, (train_idx, test_idx) in enumerate(
        splitter.split(np.zeros(len(raw_samples)), groups=groups), start=1
    ):
        fold_rows, fold_summary = _evaluate_outer_fold(
            [raw_samples[index] for index in train_idx],
            [raw_samples[index] for index in test_idx],
            ridge_grid,
            args.inner_folds,
            args.seed,
            fold,
        )
        rows.extend(fold_rows)
        folds.append(fold_summary)

    arm_names = tuple(rows[0]["arms"])
    patient_metrics: dict[str, dict[str, float]] = {}
    for arm in arm_names:
        for metric in ("average_precision", "auc", "brier"):
            values = [row["arms"][arm][metric] for row in rows]
            patient_metrics[f"{arm}.{metric}"] = {"mean_patient_value": _mean_patient_metric(rows, values)}
    patient_count_mae = _mean_patient_metric(rows, [row["count_absolute_error"] for row in rows])
    cohort_count_mae = _mean_patient_metric(rows, [row["cohort_count_absolute_error"] for row in rows])
    bootstrap = _paired_patient_bootstrap(rows, "hybrid", "matched_logistic_no_pde", "average_precision", args.seed)
    output = {
        "status": "exploratory_patient_grouped_cv",
        "independent_confirmation": False,
        "reason": "MU-Glioma-Post cases were examined in earlier Track B work; grouped CV does not create a fresh cohort.",
        "target": "newly annotated label-1/3 voxels within the forecast-origin candidate region",
        "target_limitation": "Segmentation transitions are not verified viable tumor growth or clinical recurrence; treatment effects are not modeled, and only label-4 cavity voxels are explicitly blocked (there is no CSF mask yet).",
        "cohort": {
            "n_triplets": len(raw_samples),
            "n_patients": len({sample.patient_id for sample in raw_samples}),
            "gbm_only": True,
            "patient_manifest_sha256": hashlib.sha256(args.patient_manifest.read_bytes()).hexdigest(),
            "pair_manifest_sha256": hashlib.sha256(args.pair_manifest.read_bytes()).hexdigest(),
        },
        "configuration": {
            "outer_folds": args.outer_folds,
            "inner_folds": args.inner_folds,
            "ridge_grid": list(ridge_grid),
            "seed": args.seed,
            "edge_sensitivity": diffusion.edge_sensitivity,
            "screening_length_mm": diffusion.screening_length_mm,
            "label2_log_permeability": diffusion.label2_log_permeability,
            "candidate_radius_mm": diffusion.candidate_radius_mm,
            "prior_smoothing_mm": diffusion.prior_smoothing_mm,
            "cavity_blocked": diffusion.block_cavity,
            "stability_strength": 1.0,
        },
        "primary_comparison": {
            "hybrid_minus_matched_logistic_average_precision": bootstrap,
            "hybrid_minus_hist_gradient_boosting_average_precision": _paired_patient_bootstrap(
                rows, "hybrid", "hist_gradient_boosting_no_pde", "average_precision", args.seed + 1
            ),
            "interpretation": "Exploratory patient-bootstrap interval; not evidence of external generalization or novelty.",
        },
        "count_prediction": {
            "history_shrinkage_patient_mean_absolute_error": patient_count_mae,
            "cohort_rate_only_patient_mean_absolute_error": cohort_count_mae,
        },
        "mean_patient_metrics": patient_metrics,
        "folds": folds,
        "sample_metrics": rows,
    }
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(f"refusing to overwrite {args.output}; use --overwrite intentionally")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    safe_output = _json_safe(output)
    args.output.write_text(json.dumps(safe_output, indent=2, allow_nan=False), encoding="utf-8")
    return safe_output


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument("--patient-manifest", type=Path, default=PATIENT_MANIFEST)
    parser.add_argument("--pair-manifest", type=Path, default=PAIR_MANIFEST)
    parser.add_argument("--output", type=Path, default=ROOT / "output/track_b_growth_zone_exploratory_cv.json")
    parser.add_argument("--outer-folds", type=int, default=5)
    parser.add_argument("--inner-folds", type=int, default=3)
    parser.add_argument("--ridge-grid", type=float, nargs="+", default=[0.1, 1.0, 10.0])
    parser.add_argument("--edge-sensitivity", type=float, default=0.5)
    parser.add_argument("--screening-length-mm", type=float, default=10.0)
    parser.add_argument("--label2-log-permeability", type=float, default=0.69314718056)
    parser.add_argument("--candidate-radius-mm", type=float, default=20.0)
    parser.add_argument("--prior-smoothing-mm", type=float, default=2.0)
    parser.add_argument("--open-cavity-ablation", action="store_true")
    parser.add_argument("--seed", type=int, default=20261024)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--run-exploratory-patient-cv",
        action="store_true",
        help="acknowledge that this is exploratory reuse of previously examined MU data",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    if not args.run_exploratory_patient_cv:
        raise SystemExit("Refusing to score patient data without --run-exploratory-patient-cv")
    output = run(args)
    print(f"Saved exploratory patient-grouped CV to {args.output}")
    print(f"Patients: {output['cohort']['n_patients']}; triplets: {output['cohort']['n_triplets']}")
    print(f"Hybrid - matched logistic AP (patient bootstrap): {output['primary_comparison']['hybrid_minus_matched_logistic_average_precision']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
