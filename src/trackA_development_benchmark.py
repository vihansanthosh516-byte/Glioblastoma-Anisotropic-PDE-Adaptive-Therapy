"""Run reproducible Track A implementation checks on generated inputs only."""

import argparse
import json
from pathlib import Path

import numpy as np

from .trackA_margin_evaluation import evaluate
from .trackA_regional_model import calibrate_margin_bound, fit_regional
from .trackA_tmz_model import Cultures, benchmark_tmz, matched_signature_weights


DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "output/trackA_development_benchmark/result.json"
SEED = 202611231


def _regional_fixture(rng, count):
    orientation = rng.choice([-1, 1], size=count)
    nonspatial = rng.normal(size=(count, 3)) + [0.35, 0.15, 0.45]
    spatial = rng.normal(loc=0.9, scale=0.7, size=count)
    layers = {
        "composition": nonspatial[:, 0],
        "density_depth": nonspatial[:, 1],
        "malignant_expression": nonspatial[:, 2],
        "nominal_spatial": spatial + rng.normal(scale=0.08, size=count),
        "alternative_spatial": spatial + rng.normal(scale=0.08, size=count),
        "molecule_spatial": spatial + rng.normal(scale=1.2, size=count),
    }
    features = {name: (column * orientation)[:, None] for name, column in layers.items()}
    # The generated reference is observed independently of both noisy masks.
    reference = {**features, **{
        name: (spatial * orientation)[:, None] for name in ("nominal_spatial", "alternative_spatial")
    }}
    return features, orientation, reference


def _offset(features, orientation, nominal, alternative):
    return {
        **features,
        "nominal_spatial": features["nominal_spatial"] + nominal * orientation[:, None],
        "alternative_spatial": features["alternative_spatial"] + alternative * orientation[:, None],
    }


def regional_check(rng):
    training, train_order, _ = _regional_fixture(rng, 160)
    calibration, cal_order, cal_reference = _regional_fixture(rng, 160)
    evaluation, eval_order, eval_reference = _regional_fixture(rng, 160)
    model = fit_regional(training, train_order, strength=0.1)
    calibration_reference = model.reconstruction_margins(cal_reference)[:, 0]
    truth = model.reconstruction_margins(eval_reference)[:, 0]
    conditions = {
        "clean": (0.0, 0.0, False),
        "divergent_error": (-1.8, 1.8, False),
        "shared_error": (-1.8, -1.8, False),
        "unqualified_shift": (-1.8, -1.8, True),
        "mismatched_calibration": (-1.8, -1.8, False),
    }
    candidate = "bounded_spatial"
    report = {}
    for condition, (nominal, alternative, unqualified) in conditions.items():
        matched = condition != "mismatched_calibration"
        controls = _offset(calibration, cal_order, nominal, alternative) if matched else calibration
        bound = calibrate_margin_bound(model.reconstruction_margins(controls),
                                       calibration_reference, alpha=0.05)
        measured = _offset(evaluation, eval_order, nominal, alternative)
        audit = np.full(len(eval_order), not unqualified, dtype=bool)
        margins, used = model.predict(measured, bound=bound, audit_qualified=audit)
        # Convert canonical B-A margins to the existing core/infiltrating
        # score interface for the generated patients only.
        applied = {
            name: {str(i): (0.0, float(margin * order))
                   for i, (margin, order) in enumerate(zip(values, eval_order, strict=True))}
            for name, values in margins.items()
        }
        comparisons = evaluate(applied, candidate, [arm for arm in applied if arm != candidate])
        error = np.max(np.abs(model.reconstruction_margins(measured) - truth[:, None]), axis=1)
        report[condition] = {
            "radius": bound,
            "reference_coverage": float(np.mean(error <= bound)),
            "concordance": comparisons["concordance"],
            "spatial_use": {name: float(uses.mean()) for name, uses in used.items()},
            "per_patient_concordance": comparisons["per_donor_concordance"],
            "per_patient_spatial_use": {name: uses.tolist() for name, uses in used.items()},
        }
    return report


def _cultures(rng, prefix, count, *, signal):
    sizes = rng.integers(1, 5, size=count)
    groups = np.repeat(np.arange(count), sizes)
    origins = [f"{prefix}-{i:03d}" for i in groups]
    ids = [f"{prefix}-culture-{i:03d}" for i in range(len(groups))]
    common = rng.normal(size=(count, 8))
    noise = rng.normal(size=(len(groups), 8))
    mgmt = (common[:, 7] > 0).astype(int)[groups]
    x = 0.35 * common[groups] + noise
    coefficients = np.array([0.9, -0.6, 0, 0, 0, 0, 0, 0])
    response = signal * (x @ coefficients)
    y = response + 0.22 * mgmt + rng.normal(scale=0.55, size=count)[groups]
    y += rng.normal(scale=0.35, size=len(groups))
    return Cultures(x, y, origins, ids), mgmt


def tmz_check(rng):
    fixed = np.array([1.0, -1.0, 0, 0, 0, 0, 0, 0])
    bins = ["common"] * len(fixed)
    controls = matched_signature_weights(fixed, bins, count=8, seed=SEED + 1)
    signatures = {"fixed_signature": fixed, **{f"matched_signature_{i}": x for i, x in enumerate(controls)}}
    training, train_mgmt = _cultures(rng, "train", 40, signal=1)
    evaluation, eval_mgmt = _cultures(rng, "eval", 100, signal=1)
    shifted_expression = evaluation.expression.copy()
    shifted_expression[:, 2] += 1.5 * shifted_expression[:, 0]
    shifted = Cultures(shifted_expression, evaluation.viability, evaluation.origins, evaluation.culture_ids)
    null_train, null_train_mgmt = _cultures(rng, "null_train", 40, signal=0)
    null_eval, null_eval_mgmt = _cultures(rng, "null_eval", 100, signal=0)
    conditions = {
        "signal": (training, evaluation, train_mgmt, eval_mgmt),
        "assay_shift": (training, shifted, train_mgmt, eval_mgmt),
        "null": (null_train, null_eval, null_train_mgmt, null_eval_mgmt),
    }
    report = {}
    for name, (training, evaluation, train_mgmt, eval_mgmt) in conditions.items():
        report[name] = benchmark_tmz(
            training, evaluation, penalty=0.2, within_weight=1.0,
            strata=["simulated_assay"] * len(training.origins), shuffle_seed=SEED + 2,
            mgmt=(train_mgmt, eval_mgmt), signatures=signatures,
        )
    return report


def run():
    regional_rng, tmz_rng = [np.random.default_rng(seed) for seed in np.random.SeedSequence(SEED).spawn(2)]
    regional, tmz = regional_check(regional_rng), tmz_check(tmz_rng)
    shared = regional["shared_error"]
    no_audit = regional["unqualified_shift"]
    mismatch = regional["mismatched_calibration"]
    checks = {
        "clean_spatial_is_used": regional["clean"]["spatial_use"]["bounded_spatial"] > 0,
        "shared_error_exposes_agreement": shared["concordance"]["agreement_only"] < shared["concordance"]["nonspatial"],
        "shared_error_bound_changes_use": shared["spatial_use"]["bounded_spatial"] < shared["spatial_use"]["agreement_only"],
        "unqualified_condition_falls_back": no_audit["spatial_use"]["bounded_spatial"] == 0
        and no_audit["per_patient_concordance"]["bounded_spatial"] == no_audit["per_patient_concordance"]["nonspatial"],
        "mismatched_calibration_loses_coverage": mismatch["reference_coverage"] < shared["reference_coverage"],
        "within_origin_term_changes_predictions": tmz["signal"]["within_origin"]["per_origin"]
        != tmz["signal"]["origin_mean_only"]["per_origin"],
        "sister_shuffle_changes_predictions": tmz["signal"]["within_origin"]["per_origin"]
        != tmz["signal"]["shuffled_sister_links"]["per_origin"],
    }
    return {"scope": "synthetic implementation checks only", "seed": SEED,
            "checks": checks, "regional": regional, "tmz": tmz}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"Saved synthetic implementation checks to {args.output}")
    if not all(result["checks"].values()):
        parser.exit(1, "An implementation check failed. The full output was retained.\n")


if __name__ == "__main__":
    main()
