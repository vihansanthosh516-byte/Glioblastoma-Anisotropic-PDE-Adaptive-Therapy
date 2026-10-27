"""Run the amended exploratory Track B UCSD analysis once, in its own output folder."""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
AMENDMENT = ROOT / "paper/track_b_ucsd_exploratory_amendment_2026-10-27.md"
RUNNER_PATH = ROOT / "src/166_trackb_ucsd_external_test.py"
OUTPUT_DIR = ROOT / "output/track_b_ucsd_exploratory_20261027"
MARKER_PATH = OUTPUT_DIR / "analysis_started.json"
RESULT_PATH = OUTPUT_DIR / "result.json"
ORIGINAL_MARKER = ROOT / "output/track_b_ucsd_locked_test/analysis_started.json"


def _load_runner() -> Any:
    spec = importlib.util.spec_from_file_location("trackb_ucsd_locked_runner_for_exploration", RUNNER_PATH)
    if spec is None or spec.loader is None:
        raise ImportError("could not load the frozen Track B UCSD runner")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


RUNNER = _load_runner()


def _relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT.resolve())).replace("\\", "/")


def _assert_exploratory_inputs_committed() -> tuple[str, dict[str, str]]:
    paths = (
        ".gitignore",
        _relative(AMENDMENT),
        _relative(Path(__file__)),
        "tests/test_trackb_ucsd_exploratory.py",
        *RUNNER.LOCKED_INPUT_PATHS,
    )
    status = subprocess.run(
        ["git", "status", "--porcelain", "--", *paths],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    if status:
        raise RuntimeError("exploratory protocol or scoring code is uncommitted; commit it before reading UCSD outcomes")
    ignored = subprocess.run(
        ["git", "check-ignore", "-q", _relative(RESULT_PATH)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if ignored.returncode != 0:
        raise RuntimeError("exploratory result path must be Git-ignored before UCSD outcome access")
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.strip()
    digests = {path: RUNNER._sha256(ROOT / path) for path in paths}
    return commit, digests


def _screen_origin_pairs(pairs: Sequence[Any]) -> tuple[list[tuple[int, Any]], dict[str, int], list[dict[str, Any]]]:
    """Select eligible pairs using forecast-origin MRI only; never access labels."""
    eligible: list[tuple[int, Any]] = []
    modality_failures: Counter[str] = Counter()
    failed_pairs: list[dict[str, Any]] = []
    for pair_order, pair in enumerate(pairs, start=1):
        try:
            images = RUNNER._load_ucsd_origin_images(pair.first_files)
            issues = RUNNER._origin_intensity_quality_issues(images)
        except (OSError, ValueError, RuntimeError) as exc:
            issues = ["origin_image_load_error"]
            failed_pairs.append({"pair_order": pair_order, "failure_type": type(exc).__name__})
        if issues:
            modality_failures.update(issues)
            if not any(row["pair_order"] == pair_order for row in failed_pairs):
                failed_pairs.append({"pair_order": pair_order, "modalities": sorted(set(issues))})
            continue
        eligible.append((pair_order, pair))
    return eligible, dict(sorted(modality_failures.items())), failed_pairs


def _selected_source_paths(pairs: Sequence[Any]) -> list[Path]:
    return [
        path
        for pair in pairs
        for bundle in (pair.first_files, pair.second_files)
        for path in bundle.paths.values()
    ]


def _patient_record(
    pair_order: int,
    sample: Any,
    scores: Mapping[str, np.ndarray],
    comparator: str,
) -> dict[str, Any]:
    positive = int(np.count_nonzero(sample.target == 1))
    negative = int(sample.target.size - positive)
    aucs = {name: RUNNER._auc(sample.target, score) for name, score in scores.items()}
    full_auc = aucs.get("full_hybrid")
    comparator_auc = aucs.get(comparator)
    difference = (
        float(full_auc - comparator_auc)
        if full_auc is not None and comparator_auc is not None
        else None
    )
    return {
        "pair_order": pair_order,
        "scan_gap_days": sample.interval_days,
        "scan_gap_stratum": RUNNER._interval_stratum(sample.interval_days),
        "positive_target_voxels_in_risk_set": positive,
        "negative_target_voxels_in_risk_set": negative,
        "annotated_additions_total": sample.additions_total,
        "annotated_additions_outside_risk_set": sample.additions_outside_risk,
        "annotated_additions_outside_current_support": sample.additions_outside_current_support,
        "auc_by_arm": aucs,
        "paired_auc_difference": difference,
    }


def _write_json_exclusive(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def _update_marker(**updates: Any) -> dict[str, Any]:
    marker = json.loads(MARKER_PATH.read_text(encoding="utf-8"))
    marker.update(updates)
    temporary = MARKER_PATH.with_name(MARKER_PATH.name + ".tmp")
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(marker, handle, indent=2, allow_nan=False)
        handle.write("\n")
    temporary.replace(MARKER_PATH)
    return marker


def _gap_summary(
    rows: Sequence[Mapping[str, Any]],
    eligible_pairs: Sequence[tuple[int, Any]],
) -> dict[str, Any]:
    all_strata: dict[str, list[Any]] = {key: [] for key in ("le_90_days", "91_180_days", "gt_180_days", "unknown")}
    for _, pair in eligible_pairs:
        all_strata[RUNNER._interval_stratum(pair.interval_days)].append(pair)
    row_strata: dict[str, list[Mapping[str, Any]]] = {key: [] for key in all_strata}
    for row in rows:
        row_strata[str(row["scan_gap_stratum"])].append(row)
    result: dict[str, Any] = {}
    for name in all_strata:
        paired = [row["paired_auc_difference"] for row in row_strata[name] if row["paired_auc_difference"] is not None]
        gaps = [pair.interval_days for pair in all_strata[name] if pair.interval_days is not None]
        result[name] = {
            "n_quality_eligible": len(all_strata[name]),
            "n_scored": len(row_strata[name]),
            "n_metric_evaluable": len(paired),
            "mean_paired_auc_difference": float(np.mean(paired)) if paired else None,
            "observed_gap_days": {
                "n": len(gaps),
                "median": float(np.median(gaps)) if gaps else None,
                "min": float(np.min(gaps)) if gaps else None,
                "max": float(np.max(gaps)) if gaps else None,
            },
        }
    return result


def _preflight(args: argparse.Namespace) -> tuple[
    list[Any], dict[str, Any], dict[str, int], list[dict[str, Any]], list[tuple[int, Any]]
]:
    pairs, counts = RUNNER.select_ucsd_pairs(args.ucsd_clinical, args.ucsd_root)
    if not pairs:
        raise ValueError("no GBM patient has a complete, aligned pair with all four standard modalities and masks")
    eligible, modality_counts, failed_pairs = _screen_origin_pairs(pairs)
    summary = {
        "diagnosis_eligible_patients": counts.get("gbm_patients_in_clinical_metadata", 0),
        "repeat_eligible_patients": counts.get("gbm_patients_with_repeat", 0),
        "selected_file_complete_aligned_pairs": len(pairs),
        "patients_without_file_complete_aligned_pair": counts.get("gbm_no_complete_aligned_adjacent_pair", 0),
        "origin_quality_eligible_pairs": len(eligible),
        "origin_quality_technical_failures": len(failed_pairs),
        "quality_failure_modalities": modality_counts,
        "source_pair_selection_counts": counts,
    }
    return pairs, summary, modality_counts, failed_pairs, eligible


def run_exploratory(args: argparse.Namespace) -> dict[str, Any]:
    if OUTPUT_DIR.exists():
        raise FileExistsError("exploratory output already exists; preserve it and do not overwrite")
    if not ORIGINAL_MARKER.is_file():
        raise FileNotFoundError("the original locked-test marker is missing; preserve the one-shot evidence")
    if not args.amendment.is_file() or args.amendment.resolve() != AMENDMENT.resolve():
        raise FileNotFoundError("the dated exploratory protocol amendment is required")
    if not args.patient_manifest.is_file() or not args.pair_manifest.is_file():
        raise FileNotFoundError("the frozen MU development manifests are required")
    RUNNER._validate_locked_source_paths(args)
    commit, code_digests = _assert_exploratory_inputs_committed()

    pairs, flow, _modality_counts, failed_pairs, eligible_pairs = _preflight(args)
    if not eligible_pairs:
        raise ValueError("no pair passed forecast-origin MRI quality screening")
    all_source_paths = _selected_source_paths(pairs)
    analysis_pairs = [pair for _, pair in eligible_pairs]
    ucsd_source_paths = _selected_source_paths(analysis_pairs)
    all_path_manifest = RUNNER._hash_path_size_manifest(all_source_paths, args.ucsd_root)
    path_manifest = RUNNER._hash_path_size_manifest(ucsd_source_paths, args.ucsd_root)

    parameters = RUNNER.FEATURES.OriginDiffusionParameters(
        edge_sensitivity=0.5,
        screening_length_mm=10.0,
        label2_log_permeability=0.69314718056,
        block_cavity=True,
    )
    print("loading MU development cohort", flush=True)
    development, development_source_paths = RUNNER.load_mu_development_samples(
        args.mu_data_root, args.patient_manifest, args.pair_manifest, parameters
    )
    patient_manifest_sha256 = RUNNER._sha256(args.patient_manifest)
    pair_manifest_sha256 = RUNNER._sha256(args.pair_manifest)
    development_cv, primary_comparator = RUNNER.develop_and_select_baseline(development)
    development_source_hash = RUNNER._hash_files(development_source_paths, args.mu_data_root)
    gate_scale = RUNNER._patient_gate_scale(development)
    final_models = {
        "full_hybrid": RUNNER._fit_head(
            development, scale=gate_scale, use_prior=True, use_gate=True, seed=RUNNER.FIT_SEED
        ),
        "ungated_hybrid": RUNNER._fit_head(
            development, scale=gate_scale, use_prior=True, use_gate=False, seed=RUNNER.FIT_SEED
        ),
        "mri_only": RUNNER._fit_head(
            development, scale=gate_scale, use_prior=False, use_gate=False, seed=RUNNER.FIT_SEED
        ),
    }

    output_original_sha = RUNNER._sha256(ORIGINAL_MARKER)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=False)
    marker = {
        "status": "exploratory_ucsd_outcome_access_started",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": commit,
        "amendment": _relative(AMENDMENT),
        "amendment_sha256": RUNNER._sha256(AMENDMENT),
        "scoring_code_sha256": code_digests,
        "original_locked_marker_sha256": output_original_sha,
        "original_locked_output_preserved": True,
        "cohort_flow": flow,
        "technical_failures_before_scoring": failed_pairs,
        "technical_failure_policy": "exclude origins with any nonfinite or unusable robust MRI normalization; no missing-modality fallback",
        "ucsd_all_selected_path_size_manifest_sha256": all_path_manifest,
        "ucsd_scored_path_size_manifest_sha256": path_manifest,
        "clinical_workbook_sha256": RUNNER._sha256(args.ucsd_clinical),
        "development_patient_manifest_sha256": patient_manifest_sha256,
        "development_pair_manifest_sha256": pair_manifest_sha256,
        "mu_development_source_content_manifest_sha256": development_source_hash,
        "software_versions": RUNNER._software_versions(),
        "note": "UCSD outcomes were partially exposed by the prior failed locked run; all results are exploratory.",
    }
    _write_json_exclusive(MARKER_PATH, marker)

    # Hash selected eligible source files before loading their segmentation arrays.
    ucsd_content_hash = RUNNER._hash_files(ucsd_source_paths, args.ucsd_root)
    _update_marker(ucsd_selected_content_manifest_sha256=ucsd_content_hash)

    records: list[dict[str, Any]] = []
    scoring_failures: list[dict[str, Any]] = []
    print(f"scoring {len(eligible_pairs)} quality-passing patient pairs (exploratory)", flush=True)
    for position, (pair_order, pair) in enumerate(eligible_pairs, start=1):
        try:
            current, future_labels, future_affine = RUNNER._read_pair(pair)
            sample = RUNNER._patient_sample(
                pair.patient_id, current, future_labels, future_affine, parameters, pair.interval_days
            )
            scores = RUNNER._scores(sample, scale=gate_scale, models=final_models)
            records.append(_patient_record(pair_order, sample, scores, primary_comparator))
        except Exception as exc:
            scoring_failures.append({"pair_order": pair_order, "failure_type": type(exc).__name__})
        _update_marker(
            status="exploratory_scoring_in_progress",
            pairs_processed=position,
            scoring_technical_failures=scoring_failures,
        )
        print(f"processed eligible pair {position}/{len(eligible_pairs)}", flush=True)

    if RUNNER._hash_files(ucsd_source_paths, args.ucsd_root) != ucsd_content_hash:
        raise RuntimeError("selected UCSD source files changed during exploratory scoring")
    if RUNNER._sha256(args.ucsd_clinical) != RUNNER.EXPECTED_UCSD_CLINICAL_SHA256:
        raise RuntimeError("UCSD clinical workbook changed during exploratory scoring")
    if RUNNER._sha256(args.patient_manifest) != patient_manifest_sha256:
        raise RuntimeError("MU patient manifest changed during exploratory scoring")
    if RUNNER._sha256(args.pair_manifest) != pair_manifest_sha256:
        raise RuntimeError("MU pair manifest changed during exploratory scoring")
    if RUNNER._hash_files(development_source_paths, args.mu_data_root) != development_source_hash:
        raise RuntimeError("MU development sources changed during exploratory scoring")
    if RUNNER._sha256(ORIGINAL_MARKER) != output_original_sha:
        raise RuntimeError("the original locked-test marker changed during exploratory scoring")
    if not records:
        _update_marker(status="exploratory_scoring_failed", pairs_processed=len(eligible_pairs))
        raise RuntimeError("no eligible pairs produced a score; preserve the exploratory marker and investigate")

    patient_differences = [
        {
            "pair_order": row["pair_order"],
            "scan_gap_stratum": row["scan_gap_stratum"],
            "full_hybrid_auc": row["auc_by_arm"].get("full_hybrid"),
            "comparator_auc": row["auc_by_arm"].get(primary_comparator),
            "paired_auc_difference": row["paired_auc_difference"],
        }
        for row in records
        if row["paired_auc_difference"] is not None
    ]
    differences = [row["paired_auc_difference"] for row in patient_differences]
    ci_low, ci_high = RUNNER._bootstrap_mean_ci(differences, seed=RUNNER.FIT_SEED + 100)
    auc_by_arm: dict[str, list[float]] = {"full_hybrid": []}
    auc_by_arm.update({name: [] for name in RUNNER.BASELINE_TIE_ORDER})
    for row in records:
        for name, value in row["auc_by_arm"].items():
            if value is not None:
                auc_by_arm[name].append(float(value))
    mean_auc = {name: float(np.mean(values)) if values else None for name, values in auc_by_arm.items()}
    positive_patients = sum(row["positive_target_voxels_in_risk_set"] > 0 for row in records)
    negative_patients = sum(row["negative_target_voxels_in_risk_set"] > 0 for row in records)
    scored_gaps = [row["scan_gap_days"] for row in records if row["scan_gap_days"] is not None]
    all_quality_gaps = [pair.interval_days for _, pair in eligible_pairs if pair.interval_days is not None]
    result = {
        "status": "exploratory_analysis_complete",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "interpretation": "ranking next-visit changes in source-provided NETC+ET annotations; not biological growth or clinical outcome",
        "exploratory_reason": "the original locked UCSD run accessed target labels for at least five pairs before stopping and wrote no aggregate result",
        "source": {
            "dataset": "UCSD-PTGBM v3",
            "license": "CC BY 4.0 per the current TCIA collection page",
            "dataset_doi": "10.7937/fwv2-dt74",
            "dataset_article": "Gagnon et al., Scientific Data (2025), doi:10.1038/s41597-025-06499-z",
            "clinical_workbook_sha256": RUNNER._sha256(args.ucsd_clinical),
            "development_patient_manifest_sha256": patient_manifest_sha256,
            "development_pair_manifest_sha256": pair_manifest_sha256,
            "all_selected_source_path_size_manifest_sha256": all_path_manifest,
            "scored_source_path_size_manifest_sha256": path_manifest,
            "selected_source_content_manifest_sha256": ucsd_content_hash,
            "mu_development_source_content_manifest_sha256": development_source_hash,
            "original_locked_marker_sha256": output_original_sha,
            "software_versions": RUNNER._software_versions(),
        },
        "development": {
            "dataset": "MU-Glioma-Post",
            "n_patients": len(development),
            "patient_grouped_cv_auc_for_baseline_selection": development_cv,
            "primary_comparator_selected_without_UCSD_outcomes": primary_comparator,
            "gate_scale_from_development_only": gate_scale,
        },
        "cohort": {
            **flow,
            "scored_pairs": len(records),
            "scoring_technical_failures": scoring_failures,
            "metric_evaluable_patients": len(differences),
            "patients_with_positive_target_voxels": int(positive_patients),
            "patients_with_negative_target_voxels": int(negative_patients),
            "patients_with_zero_in_risk_additions": sum(
                row["positive_target_voxels_in_risk_set"] == 0 for row in records
            ),
            "positive_target_voxels_in_risk_set": int(
                sum(row["positive_target_voxels_in_risk_set"] for row in records)
            ),
            "negative_target_voxels_in_risk_set": int(
                sum(row["negative_target_voxels_in_risk_set"] for row in records)
            ),
            "annotated_additions_total": int(sum(row["annotated_additions_total"] for row in records)),
            "annotated_additions_outside_risk_set": int(
                sum(row["annotated_additions_outside_risk_set"] for row in records)
            ),
            "annotated_additions_outside_current_support": int(
                sum(row["annotated_additions_outside_current_support"] for row in records)
            ),
            "quality_eligible_scan_gap_days": {
                "n": len(all_quality_gaps),
                "median": float(np.median(all_quality_gaps)) if all_quality_gaps else None,
                "min": float(np.min(all_quality_gaps)) if all_quality_gaps else None,
                "max": float(np.max(all_quality_gaps)) if all_quality_gaps else None,
            },
            "scored_scan_gap_days": {
                "n": len(scored_gaps),
                "median": float(np.median(scored_gaps)) if scored_gaps else None,
                "min": float(np.min(scored_gaps)) if scored_gaps else None,
                "max": float(np.max(scored_gaps)) if scored_gaps else None,
            },
            "scan_gap_strata": _gap_summary(records, eligible_pairs),
        },
        "primary_comparison": {
            "contrast": f"full_hybrid minus MU-selected {primary_comparator}",
            "n_evaluable_patients": len(differences),
            "mean_paired_patient_auc_difference": float(np.mean(differences)) if differences else None,
            "paired_patient_bootstrap_95_percentile_ci": [ci_low, ci_high] if differences else None,
            "bootstrap_draws": RUNNER.N_BOOTSTRAPS,
            "bootstrap_seed": RUNNER.FIT_SEED + 100,
            "former_locked_protocol_margin_auc": RUNNER.MIN_EFFECT_AUC,
            "confirmatory_success_criterion_claimed": False,
        },
        "patient_macro_auc": {
            name: {"n_evaluable_patients": len(values), "mean_auc": mean_auc[name]}
            for name, values in auc_by_arm.items()
        },
        "per_patient_paired_auc_differences": patient_differences,
        "per_patient_summary": records,
        "technical_failures_before_scoring": failed_pairs,
        "analysis": {
            "patient_identifiers_written": False,
            "voxel_predictions_written": False,
            "missing_modality_fallback_used": False,
            "scan_gap_used_as_model_input": False,
            "one_shot_marker": _relative(MARKER_PATH),
            "result_file": _relative(RESULT_PATH),
            "model_code": "unchanged src/166_trackb_ucsd_external_test.py and its frozen dependencies",
        },
        "limitations": [
            "This is exploratory because UCSD target labels were partially accessed in the failed locked run.",
            "The MRI-only quality exclusion can select a nonrepresentative subset of eligible patients.",
            "The endpoint is a change in supplied annotations; treatment, registration, and contouring can affect it.",
            "Scan gaps vary and are not model inputs; results do not estimate a fixed-horizon forecast.",
            "The boundary gate is a fixed sensitivity heuristic, not calibrated uncertainty.",
            "Published-method comparisons are contextual; no head-to-head superiority claim is made.",
        ],
    }
    _write_json_exclusive(RESULT_PATH, result)
    _update_marker(
        status="exploratory_analysis_complete",
        result_sha256=RUNNER._sha256(RESULT_PATH),
        scored_pairs=len(records),
        metric_evaluable_patients=len(differences),
        scoring_technical_failures=len(scoring_failures),
    )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ucsd-root", type=Path, default=RUNNER.DEFAULT_UCSD_ROOT)
    parser.add_argument("--ucsd-clinical", type=Path, default=RUNNER.DEFAULT_UCSD_CLINICAL)
    parser.add_argument("--mu-data-root", type=Path, default=RUNNER.DEFAULT_MU_ROOT)
    parser.add_argument("--patient-manifest", type=Path, default=RUNNER.DEFAULT_PATIENT_MANIFEST)
    parser.add_argument("--pair-manifest", type=Path, default=RUNNER.DEFAULT_PAIR_MANIFEST)
    parser.add_argument("--protocol", type=Path, default=RUNNER.DEFAULT_PROTOCOL)
    parser.add_argument("--amendment", type=Path, default=AMENDMENT)
    parser.add_argument("--preflight-only", action="store_true", help="screen origin MRI only; do not read any segmentation voxels")
    parser.add_argument("--score-exploratory", action="store_true", help="run the dated exploratory analysis once")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.preflight_only == args.score_exploratory:
        parser.error("choose exactly one: --preflight-only or --score-exploratory")
    if args.preflight_only:
        RUNNER._validate_locked_source_paths(args)
        _pairs, summary, _modalities, failures, _eligible = _preflight(args)
        print(json.dumps({**summary, "technical_failure_records": failures, "ucsd_segmentation_voxels_read": False}, indent=2))
        return 0
    result = run_exploratory(args)
    print(json.dumps({
        "status": result["status"],
        "scored_pairs": result["cohort"]["scored_pairs"],
        "metric_evaluable_patients": result["cohort"]["metric_evaluable_patients"],
        "mean_paired_auc_difference": result["primary_comparison"]["mean_paired_patient_auc_difference"],
        "result_path": result["analysis"]["result_file"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
