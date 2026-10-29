"""Training-only volume nuisance fitting for MU three-visit development folds.

This module reads only the three visits in supplied outer-training records. It
does not load files or inspect held-out patients.
"""
from __future__ import annotations

import numpy as np

from src.fk_patient_cv import PatientSeries
from src.fk_volume_model import (
    fit_log_volume_variance,
    select_shrinkage_time_constant,
)


def fit_mu_training_volume_nuisance(
    training_records, *, tau_candidates_days=(30, 90, 180, 365)
) -> dict:
    """Fit population slope, shrinkage time, and horizon variance on training patients.

    Each record must contain exactly three visits. Patients with any empty or
    nonpositive mask volume are excluded from nuisance fitting and counted.
    The caller keeps such patients in outer CV scoring and reports failures.
    """
    records = tuple(training_records)
    if not records:
        raise ValueError("training_records must contain exactly-three-visit PatientSeries records")

    seen_ids = set()
    valid = []
    invalid_volume_count = 0
    for record in records:
        if not isinstance(record, PatientSeries):
            raise ValueError("training_records must contain PatientSeries values")
        patient_id = str(record.patient_id).strip()
        if not patient_id or patient_id in seen_ids:
            raise ValueError("training patient IDs must be nonempty and unique")
        seen_ids.add(patient_id)
        days = np.asarray(record.days, dtype=float)
        masks = tuple(record.masks)
        if (days.shape != (3,) or len(masks) != 3 or not np.isfinite(days).all()
                or np.any(np.diff(days) <= 0)):
            raise ValueError("each training record must have exactly three ordered visits")
        volumes = []
        for mask in masks:
            array = np.asarray(mask)
            if array.ndim != 3 or not (array.dtype == np.bool_ or np.isin(array, (0, 1)).all()):
                raise ValueError("training masks must be binary 3D arrays")
            volumes.append(int(np.count_nonzero(array)))
        if any(volume <= 0 for volume in volumes):
            invalid_volume_count += 1
            continue
        valid.append((patient_id, days, np.asarray(volumes, dtype=float)))

    if len(valid) < 3:
        raise ValueError("need at least three positive-volume training patients")

    patient_slopes = {}
    for patient_id, days, volumes in valid:
        patient_slopes[patient_id] = float(
            (np.log(volumes[1]) - np.log(volumes[0])) / (days[1] - days[0])
        )
    population_slope = float(np.median(list(patient_slopes.values())))

    validation_cases = []
    for patient_id, days, volumes in valid:
        other_slopes = [slope for other_id, slope in patient_slopes.items() if other_id != patient_id]
        loo_population = float(np.median(other_slopes))
        first_gap = float(days[1] - days[0])
        horizon = float(days[2] - days[1])
        validation_cases.append({
            "patient_id": patient_id,
            "patient_slope": patient_slopes[patient_id],
            "population_slope": loo_population,
            "observation_gap_days": first_gap,
            "horizon_days": horizon,
            "last_volume": float(volumes[1]),
            "target_volume": float(volumes[2]),
        })

    selected_tau = select_shrinkage_time_constant(
        validation_cases, candidates_days=tau_candidates_days
    )
    residuals, horizons, residual_patient_ids = [], [], []
    tau = selected_tau["time_constant_days"]
    from src.fk_volume_model import shrink_log_growth_slope

    for case in validation_cases:
        shrunk_slope = shrink_log_growth_slope(
            case["patient_slope"], case["population_slope"],
            observation_gap_days=case["observation_gap_days"], time_constant_days=tau,
        )
        predicted = np.log(case["last_volume"]) + shrunk_slope * case["horizon_days"]
        residuals.append(float(np.log(case["target_volume"]) - predicted))
        horizons.append(case["horizon_days"])
        residual_patient_ids.append(case["patient_id"])

    if len(np.unique(horizons)) < 2:
        raise ValueError("training records must provide at least two effective forecast horizons")
    variance_fit = fit_log_volume_variance(
        horizons, residuals, patient_ids=residual_patient_ids
    )
    return {
        "population_log_volume_slope_per_day": population_slope,
        "shrinkage_time_constant_days": tau,
        "shrinkage_selection_errors": selected_tau["candidate_mse_log_volume"],
        "volume_variance_fit": variance_fit,
        "audit_counts": {
            "n_training_patients": len(records),
            "n_positive_volume_patients": len(valid),
            "n_excluded_nonpositive_volume_patients": invalid_volume_count,
            "n_validation_forecasts": len(validation_cases),
            "n_effective_horizons_days": len(np.unique(horizons)),
        },
        "training_patient_ids": tuple(sorted(seen_ids)),
        "scope": "outer-training patients only; exactly three visits per patient",
    }
