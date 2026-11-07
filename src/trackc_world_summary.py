"""World-level summaries for paired simulated policy differences.

This development helper does not define an adopted scoring protocol.
"""
from __future__ import annotations

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
    differences = (candidate - comparator).mean(axis=1)
    rng = np.random.default_rng(seed)
    boot = differences[rng.integers(0, n_worlds, size=(resamples, n_worlds))].mean(axis=1)
    return {
        "n_worlds": n_worlds,
        "n_patients_per_world": n_patients,
        "world_mean_differences": differences,
        "equal_world_mean_difference": float(differences.mean()),
        "bootstrap_percentile_95": np.quantile(boot, [0.025, 0.975]),
        "bootstrap_unit": "world",
        "resamples": int(resamples),
        "seed": seed,
    }


def summarize_benchmark_worlds(records, *, resamples, seed, contrast="loo"):
    """Summarize a declared LOO or finite-set minimax development contrast.

    Matching recorded settings and distinct fingerprints are necessary checks,
    not proof that independently generated worlds or patient pairings are valid.
    """
    records = tuple(records)
    if contrast not in ("loo", "minimax"):
        raise ValueError("contrast must be loo or minimax")
    if len(records) < 2:
        raise ValueError("at least two benchmark worlds are required")
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
        if contrast == "loo":
            recorded_difference = np.asarray(record["primary_paired_contrast"]["robust_utility_difference"], dtype=float)
            if record["primary_paired_contrast"]["contrast"] != "loo_dependence_aware_mpc minus loo_conditional_independence_mpc":
                raise ValueError("unexpected primary contrast")
            if recorded_difference.shape != candidate.shape or not np.array_equal(recorded_difference, candidate - comparator):
                raise ValueError("recorded primary contrast does not match paired outcomes")
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
    result["role"] = "synthetic development; not an adopted confirmatory analysis"
    result["weight_source"] = "evaluation cohort peers" if contrast == "loo" else "fixed finite model sets"
    result["checks_limit"] = "Distinct recorded seeds and hashes do not prove independence or pairing."
    return result
