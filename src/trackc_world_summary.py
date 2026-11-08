"""World-level summaries for paired simulated policy differences.

This development helper does not define an adopted scoring protocol.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np


def paired_world_summary(candidate, comparator, world_ids, *, resamples, seed):
    """Average patients within worlds, then bootstrap independent worlds equally.

    Inputs are matching two-dimensional world-by-patient arrays. Patient positions
    must represent paired simulated patients. World identifiers must be unique.
    """
    candidate = np.asarray(candidate, dtype=float)
    comparator = np.asarray(comparator, dtype=float)
    if candidate.ndim != 2 or candidate.shape != comparator.shape:
        raise ValueError("paired outcomes must be matching world-by-patient arrays")
    n_worlds, n_patients = candidate.shape
    if n_worlds < 2 or n_patients < 1:
        raise ValueError("at least two worlds and one patient per world are required")
    if not np.isfinite(candidate).all() or not np.isfinite(comparator).all():
        raise ValueError("paired outcomes must be finite")
    ids = tuple(world_ids)
    if len(ids) != n_worlds or any(not isinstance(x, str) or not x.strip() for x in ids) or len(set(ids)) != n_worlds:
        raise ValueError("world identifiers must be unique nonempty strings")
    if isinstance(resamples, (bool, np.bool_)) or not isinstance(resamples, (int, np.integer)) or resamples < 2:
        raise ValueError("resamples must be an integer of at least 2")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    differences = (candidate - comparator).mean(axis=1)
    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    boot = differences[rng.integers(0, n_worlds, size=(resamples, n_worlds))].mean(axis=1)
    return {
        "n_worlds": n_worlds,
        "n_patients_per_world": n_patients,
        "world_mean_differences": differences,
        "equal_world_mean_difference": float(differences.mean()),
        "bootstrap_percentile_95": np.quantile(boot, [0.025, 0.975]),
        "bootstrap_unit": "world",
        "resamples": int(resamples),
        "seed": int(seed),
        "rng": "numpy.random.Generator(PCG64DXSM)",
    }


def summarize_benchmark_worlds(records, *, resamples, seed, contrast=None):
    """Summarize a declared LOO or finite-set minimax development contrast.

    Matching recorded settings and distinct fingerprints are necessary checks,
    not proof that independently generated worlds or patient pairings are valid.
    When ``contrast`` is omitted, use the contrast explicitly marked primary in
    every record; legacy records without that field default to the LOO contrast.
    """
    records = tuple(records)
    if len(records) < 2:
        raise ValueError("at least two benchmark worlds are required")
    if contrast is None:
        primary_names = {record.get("primary_contrast_name", "loo") for record in records}
        if len(primary_names) != 1:
            raise ValueError("benchmark worlds have inconsistent primary contrast names")
        contrast = primary_names.pop()
    if contrast not in ("loo", "minimax"):
        raise ValueError("contrast must be loo or minimax")
    reference = records[0]
    candidates, comparators, world_ids = [], [], []
    seen_fingerprints = set()
    seen_seed_keys = set()
    seen_world_hashes = set()
    reference_family = None
    for record in records:
        if record["analysis_type"] != "synthetic_virtual_world_development":
            raise ValueError("all records must be synthetic development worlds")
        if contrast == "minimax":
            provenance = record.get("world_provenance")
            if not isinstance(provenance, dict):
                raise ValueError("minimax summary requires world provenance")
            if "scan_loss" not in provenance:
                raise ValueError("minimax summary requires scan-loss provenance")
            if provenance["scan_loss"] != reference.get("world_provenance", {}).get("scan_loss"):
                raise ValueError("world scan-loss settings differ")
            family = provenance.get("family")
            world_hash = provenance.get("world_parameters_sha256")
            world_seed = provenance.get("world_generation_seed")
            if family not in ("single", "mixture"):
                raise ValueError("invalid world family")
            if reference_family is None:
                reference_family = family
            elif family != reference_family:
                raise ValueError("world families differ")
            if not isinstance(world_hash, str) or len(world_hash) != 64 or any(c not in "0123456789abcdef" for c in world_hash):
                raise ValueError("invalid world parameter fingerprint")
            if world_hash in seen_world_hashes:
                raise ValueError("world parameter fingerprint reused")
            seen_world_hashes.add(world_hash)
            if isinstance(world_seed, (bool, np.bool_)) or not isinstance(world_seed, (int, np.integer)) or world_seed < 0:
                raise ValueError("invalid world generation seed")
        for key in ("settings", "inference_model_parameters"):
            if record[key] != reference[key]:
                raise ValueError(f"world {key} differs")
        for key in ("generator", "parameters", "f0_source", "mri_template_sha256"):
            if record["signal_generation"][key] != reference["signal_generation"][key]:
                raise ValueError(f"world signal-generation {key} differs")
        if record["signal_generation"]["f0_source"] != "evaluation_scenarios['f0']":
            raise ValueError("unexpected evaluation signal source")
        calibration_signal = record["calibration_signal_generation"]
        if calibration_signal["generator"] != record["signal_generation"]["generator"] or calibration_signal["parameters"] != record["signal_generation"]["parameters"]:
            raise ValueError("calibration and evaluation signal generation differ")
        if calibration_signal["f0_source"] != "calibration_scenarios['f0']":
            raise ValueError("unexpected calibration signal source")
        for key in ("selection_scenarios_sha256", "calibration_scenarios_sha256", "evaluation_scenarios_sha256"):
            fingerprint = record["input_fingerprints"][key]
            if not isinstance(fingerprint, str) or len(fingerprint) != 64 or any(c not in "0123456789abcdef" for c in fingerprint):
                raise ValueError("invalid scenario fingerprint")
            if fingerprint in seen_fingerprints:
                raise ValueError("scenario fingerprint reused")
            seen_fingerprints.add(fingerprint)
        seed_keys = [tuple(value) for value in record["seeds"].values()]
        seed_keys.append((record["signal_generation"]["seed"],))
        if contrast == "minimax":
            seed_keys.append((int(world_seed),))
        for seed_key in seed_keys:
            if not seed_key or seed_key in seen_seed_keys:
                raise ValueError("seed key reused or empty")
            seen_seed_keys.add(seed_key)
        arms = record["arm_outcomes"]
        arm_names = ("loo_dependence_aware_mpc", "loo_conditional_independence_mpc") if contrast == "loo" else ("finite_set_minimax_mpc", "finite_set_d0_minimax_mpc")
        if any(name not in arms for name in arm_names):
            raise ValueError(f"required {contrast} arms are missing")
        candidate = np.asarray(arms[arm_names[0]]["robust_utility"], dtype=float)
        comparator = np.asarray(arms[arm_names[1]]["robust_utility"], dtype=float)
        if candidate.ndim != 1 or candidate.shape != comparator.shape or candidate.size != record["n_evaluation_patients"]:
            raise ValueError("world outcome arrays do not match the evaluation cohort")
        expected_labels = {
            "loo": "loo_dependence_aware_mpc minus loo_conditional_independence_mpc",
            "minimax": "finite_set_minimax_mpc minus finite_set_d0_minimax_mpc",
        }
        paired = record.get("paired_contrasts", {}).get(contrast)
        if paired is None and record.get("primary_contrast_name") == contrast:
            paired = record.get("primary_paired_contrast")
        if paired is None and contrast == "loo":
            # Read legacy records that stored only the LOO primary contrast.
            paired = record.get("primary_paired_contrast")
        if paired is not None:
            if paired.get("contrast") != expected_labels[contrast]:
                raise ValueError(f"unexpected {contrast} contrast")
            if paired.get("candidate_arm") not in (None, arm_names[0]) or paired.get("control_arm") not in (None, arm_names[1]):
                raise ValueError(f"unexpected {contrast} contrast arms")
            recorded_difference = np.asarray(paired["robust_utility_difference"], dtype=float)
            if recorded_difference.shape != candidate.shape or not np.array_equal(recorded_difference, candidate - comparator):
                raise ValueError(f"recorded {contrast} paired difference does not match arm outcomes")
        candidates.append(candidate)
        comparators.append(comparator)
        world_ids.append(record["world_id"])
    try:
        candidate_matrix = np.stack(candidates)
        comparator_matrix = np.stack(comparators)
    except ValueError as exc:
        raise ValueError("worlds must have equal evaluation cohort sizes") from exc
    result = paired_world_summary(candidate_matrix, comparator_matrix, world_ids, resamples=resamples, seed=seed)
    result["contrast"] = contrast
    result["primary_contrast_name"] = contrast
    result["role"] = "synthetic development; not an adopted confirmatory analysis"
    result["weight_source"] = "evaluation cohort peers" if contrast == "loo" else "fixed finite model sets"
    result["checks_limit"] = "Distinct recorded seeds and hashes do not prove independence or pairing."
    return result


def simultaneous_max_t_world_intervals(
    world_differences_by_family,
    *,
    resamples,
    seed,
    confidence_level=0.95,
    precision_target_half_width=14.0,
):
    """Build familywise intervals from equal-weight, paired world means.

    Each input vector contains one candidate-minus-control patient-macro mean per
    independent world. Candidate and control must use the same patients within a
    world. Resampling occurs within response family at the world level. The max-T
    statistic uses each observed family standard error as its studentizer; the
    reported symmetric intervals use the higher empirical quantile of the max
    absolute centered bootstrap statistics.
    """
    if not isinstance(world_differences_by_family, Mapping) or len(world_differences_by_family) < 2:
        raise ValueError("at least two response families are required")
    if isinstance(resamples, (bool, np.bool_)) or not isinstance(resamples, (int, np.integer)) or resamples < 2:
        raise ValueError("resamples must be an integer of at least 2")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    try:
        confidence_level = float(confidence_level)
        precision_target_half_width = float(precision_target_half_width)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("confidence and precision settings must be numeric") from exc
    if not np.isfinite(confidence_level) or not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must lie strictly between 0 and 1")
    if not np.isfinite(precision_target_half_width) or precision_target_half_width <= 0.0:
        raise ValueError("precision_target_half_width must be finite and positive")

    family_names = tuple(world_differences_by_family)
    if any(not isinstance(family, str) or not family.strip() for family in family_names):
        raise ValueError("response family names must be nonempty strings")
    expected_families = {"single-response", "mixture-response"}
    if set(family_names) != expected_families:
        raise ValueError("families must be exactly single-response and mixture-response")

    family_values = {}
    family_statistics = []
    family_estimates = {}
    family_standard_errors = {}
    for index, family in enumerate(sorted(family_names)):
        values = np.asarray(world_differences_by_family[family], dtype=float)
        if values.ndim != 1 or values.size < 2 or not np.isfinite(values).all():
            raise ValueError(f"{family} must contain at least two finite world differences")
        with np.errstate(over="ignore", invalid="ignore"):
            estimate = float(values.mean())
            standard_error = float(values.std(ddof=1) / np.sqrt(values.size))
        if not np.isfinite(estimate) or not np.isfinite(standard_error) or standard_error <= 0.0:
            raise ValueError(f"{family} world differences must have nonzero variance")
        rng = np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence([int(seed), index])))
        indices = rng.integers(0, values.size, size=(int(resamples), values.size))
        bootstrap_means = values[indices].mean(axis=1)
        family_statistics.append(np.abs((bootstrap_means - estimate) / standard_error))
        family_values[family] = values.copy()
        family_estimates[family] = estimate
        family_standard_errors[family] = standard_error

    max_statistics = np.max(np.stack(family_statistics, axis=1), axis=1)
    critical_value = float(np.quantile(max_statistics, confidence_level, method="higher"))
    intervals = {}
    for family, values in family_values.items():
        estimate = family_estimates[family]
        standard_error = family_standard_errors[family]
        half_width = critical_value * standard_error
        intervals[family] = {
            "n_worlds": int(values.size),
            "world_mean_differences": values.tolist(),
            "equal_world_mean_difference": estimate,
            "standard_error": standard_error,
            "adjusted_lower": estimate - half_width,
            "adjusted_upper": estimate + half_width,
            "adjusted_half_width": half_width,
            "precision_pass": bool(half_width <= precision_target_half_width),
        }

    all_positive = all(row["adjusted_lower"] > 0.0 for row in intervals.values())
    precision_pass = all(row["precision_pass"] for row in intervals.values())
    if all_positive and precision_pass:
        decision = "pass_if_adopted"
    elif not precision_pass:
        decision = "inconclusive_precision_target_missed"
    else:
        decision = "does_not_pass_positive_rule"
    return {
        "method": "world-resampled observed-SE max-T bootstrap",
        "resampling_unit": "one paired candidate-control patient-mean difference per independent world; worlds resampled separately within each response family",
        "dependence_assumption": "World samples are independent within and across response families; fingerprints alone cannot establish this.",
        "critical_value_rule": "higher empirical quantile of the max absolute centered bootstrap statistic",
        "interval_rule": "equal-world mean plus or minus max-T critical value times observed standard error",
        "confidence_level": confidence_level,
        "resamples": int(resamples),
        "seed": int(seed),
        "protocol_status": "proposal_not_frozen",
        "rng": "numpy.random.Generator(PCG64DXSM) with family-index SeedSequence streams",
        "precision_target_half_width_unit": "simulated days",
        "max_t_critical_value": critical_value,
        "precision_target_half_width": precision_target_half_width,
        "all_adjusted_lower_bounds_above_zero": bool(all_positive),
        "precision_pass": bool(precision_pass),
        "decision": decision,
        "intervals_by_family": intervals,
        "interpretation_limit": "Synthetic response-family inference only; intervals do not validate the simulator or imply patient benefit.",
    }
