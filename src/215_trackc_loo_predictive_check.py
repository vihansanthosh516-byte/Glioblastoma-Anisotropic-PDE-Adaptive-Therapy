"""Frozen synthetic LOO paired-signal predictive-score diagnostic.

This check evaluates only synthetic baseline-signal prediction. It does not run a
treatment policy, read patient data, or make a clinical claim. See the frozen
protocol in paper/trackc_loo_predictive_check_protocol_2026-11-06.md.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

from trackc_integrated_signal_model import SignalModelParameters


ROOT = Path(__file__).resolve().parent.parent
PROTOCOL_PATH = ROOT / "paper" / "trackc_loo_predictive_check_protocol_2026-11-06.md"
ROOT_SEEDS = tuple(range(20261905, 20261930))
SAMPLE_SIZES = (18, 28, 100, 500)
N_F0_NODES = 256
NOMINAL_PARAMETERS = SignalModelParameters(0.8, 0.5, 0.4, 0.6)
PRIMARY_TRUTHS = (
    (0.7, 0.5, 0.5, -0.5),
    (0.7, 0.5, 0.5, 0.5),
    (0.9, 0.5, 0.5, 0.5),
)
ADVANCEMENT_SAMPLE_SIZES = (18, 28)
FAMILYWISE_ONE_SIDED_LEVEL = 0.9875


def _load_sibling(name: str, filename: str):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {filename}: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_RECOVERY = _load_sibling("trackc_recovery_for_loo_predictive_check", "213_trackc_finite_sample_signal_recovery.py")
_LOO = _load_sibling("trackc_loo_for_predictive_check", "214_trackc_loo_model_average.py")
_SEED_REGISTRY = _RECOVERY._SEED_REGISTRY
F0_RANGE_NAMES = _RECOVERY.F0_RANGE_NAMES
F0_RANGES = _RECOVERY.F0_RANGES
TRUTH_ANCHORS = _RECOVERY.TRUTH_ANCHORS
MODEL_GRID = tuple(_LOO.DEFAULT_MODEL_GRID)


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _file_hashes() -> dict[str, str]:
    paths = (
        Path(__file__),
        Path(__file__).with_name("214_trackc_loo_model_average.py"),
        Path(__file__).with_name("213_trackc_finite_sample_signal_recovery.py"),
        Path(__file__).with_name("212_trackc_signal_fisher_grid.py"),
        Path(__file__).with_name("211_trackc_signal_fisher_identifiability.py"),
        Path(__file__).with_name("210_trackc_dependence_identifiability.py"),
        Path(__file__).with_name("208_trackc_robust_belief_mpc.py"),
        Path(__file__).with_name("trackc_integrated_signal_model.py"),
        Path(__file__).with_name("trackc_seed_registry.py"),
        PROTOCOL_PATH,
    )
    return {str(path.relative_to(ROOT)): _hash(path) for path in paths}


def _fixed_model_score(observations, parameters, f0_bounds, n_f0_nodes):
    values = _LOO.marginal_signal_log_densities(
        observations,
        model_parameters=(parameters,),
        f0_bounds=f0_bounds,
        n_f0_nodes=n_f0_nodes,
    )
    return values[:, 0]


def _cohort_scores(observations, truth, f0_bounds, n_f0_nodes):
    predictive = _LOO.loo_signal_predictive_scores(
        observations,
        model_parameters=MODEL_GRID,
        f0_bounds=f0_bounds,
        n_f0_nodes=n_f0_nodes,
    )
    independence = SignalModelParameters(truth[0], truth[1], truth[2], 0.0)
    oracle = SignalModelParameters(*truth)
    individual = {
        "loo_model_average": predictive["loo_log_predictive_density"],
        "uniform_grid_mixture": predictive["uniform_grid_log_predictive_density"],
        "nominal_fixed": _fixed_model_score(
            observations, NOMINAL_PARAMETERS, f0_bounds, n_f0_nodes
        ),
        "matched_conditional_independence": _fixed_model_score(
            observations, independence, f0_bounds, n_f0_nodes
        ),
        "generating_parameter_oracle": _fixed_model_score(
            observations, oracle, f0_bounds, n_f0_nodes
        ),
    }
    n_patients = len(np.asarray(observations.molecular_report))
    if any(np.asarray(values).shape != (n_patients,) for values in individual.values()):
        raise ArithmeticError("a predictive-score arm returned the wrong patient count")
    if any(not np.isfinite(values).all() for values in individual.values()):
        raise ArithmeticError("a predictive-score arm returned a non-finite score")
    return {name: float(np.mean(values)) for name, values in individual.items()}


def _is_primary_truth(truth) -> bool:
    return tuple(map(float, truth)) in PRIMARY_TRUTHS


def _mean_and_lower_bound(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size < 2 or not np.isfinite(values).all():
        return {"n_root_replicates": int(values.size), "mean": None, "sd": None, "one_sided_lower_bound": None}
    critical = float(student_t.ppf(FAMILYWISE_ONE_SIDED_LEVEL, df=values.size - 1))
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=1))
    return {
        "n_root_replicates": int(values.size),
        "mean": mean,
        "sd": sd,
        "one_sided_level": FAMILYWISE_ONE_SIDED_LEVEL,
        "one_sided_lower_bound": float(mean - critical * sd / np.sqrt(values.size)),
    }


def _primary_gate(rows, sample_sizes, seed_roots, *, protocol_run):
    contrasts = {
        "vs_uniform_grid_mixture": ("loo_model_average", "uniform_grid_mixture"),
        "vs_matched_conditional_independence": (
            "loo_model_average", "matched_conditional_independence"
        ),
    }
    results = {}
    for n in ADVANCEMENT_SAMPLE_SIZES:
        if n not in sample_sizes:
            continue
        for contrast_name, (candidate, comparator) in contrasts.items():
            root_means = []
            for root in seed_roots:
                root_rows = [
                    row for row in rows
                    if row["n"] == n and row["seed_root"] == root
                    and tuple(row["truth_parameters"]) in PRIMARY_TRUTHS
                ]
                expected = len(F0_RANGE_NAMES) * len(PRIMARY_TRUTHS)
                if len(root_rows) != expected:
                    raise ArithmeticError(
                        f"expected {expected} primary rows for root {root}, n={n}; found {len(root_rows)}"
                    )
                root_means.append(
                    float(np.mean([
                        row["scores"][candidate] - row["scores"][comparator]
                        for row in root_rows
                    ]))
                )
            results[f"n={n}|{contrast_name}"] = _mean_and_lower_bound(root_means)
    required = {
        f"n={n}|{contrast}"
        for n in ADVANCEMENT_SAMPLE_SIZES
        for contrast in contrasts
    }
    present = required.issubset(results)
    passed = protocol_run and present and all(
        results[key]["one_sided_lower_bound"] is not None
        and results[key]["one_sided_lower_bound"] > 0.0
        for key in required
    )
    return {
        "status": (
            "pass" if passed else
            "not_protocol_run" if not protocol_run else
            "not_passed" if present else
            "incomplete_design"
        ),
        "rule": "all four one-sided 98.75% t lower bounds at n=18 and n=28 must be > 0",
        "familywise_one_sided_level": FAMILYWISE_ONE_SIDED_LEVEL,
        "comparisons": results,
    }


def build_diagnostic(
    *,
    seed_roots=ROOT_SEEDS,
    sample_sizes=SAMPLE_SIZES,
    n_f0_nodes=N_F0_NODES,
    check_registered_seeds=True,
):
    """Build the frozen protocol result; reduced arguments are for smoke tests only."""
    pre_run_code_hashes = _file_hashes()
    pre_run_protocol_hash = _hash(PROTOCOL_PATH)
    seed_roots = tuple(int(value) for value in seed_roots)
    sample_sizes = tuple(int(value) for value in sample_sizes)
    if not seed_roots or len(set(seed_roots)) != len(seed_roots) or any(value < 0 for value in seed_roots):
        raise ValueError("seed_roots must be unique nonnegative integers")
    if not sample_sizes or len(set(sample_sizes)) != len(sample_sizes) or any(value < 1 for value in sample_sizes):
        raise ValueError("sample_sizes must be unique positive integers")
    if isinstance(n_f0_nodes, bool) or int(n_f0_nodes) != n_f0_nodes or n_f0_nodes < 2:
        raise ValueError("n_f0_nodes must be an integer of at least 2")
    if check_registered_seeds:
        problems = _SEED_REGISTRY.check_free(
            scenario=(min(seed_roots), max(seed_roots)), purpose="stage2_dev"
        )
        if problems:
            raise RuntimeError(f"registered development seed check failed: {problems}")

    rows = []
    for range_index, range_name in enumerate(F0_RANGE_NAMES):
        f0_bounds = F0_RANGES[range_name]
        for truth_index, truth in enumerate(TRUTH_ANCHORS):
            truth_label = f"q{truth[0]}_a{truth[1]}_b{truth[2]}_d{truth[3]}"
            for size_index, n in enumerate(sample_sizes):
                for seed_root in seed_roots:
                    observations = _RECOVERY._seeded_dataset(
                        seed_root, range_index, truth_index, size_index, truth, f0_bounds, n
                    )
                    scores = _cohort_scores(observations, truth, f0_bounds, int(n_f0_nodes))
                    rows.append({
                        "f0_range": range_name,
                        "f0_bounds": list(map(float, f0_bounds)),
                        "truth_label": truth_label,
                        "truth_parameters": list(map(float, truth)),
                        "primary_informative_truth": _is_primary_truth(truth),
                        "n": int(n),
                        "seed_root": int(seed_root),
                        "scores": scores,
                        "delta_loo_vs_uniform": scores["loo_model_average"] - scores["uniform_grid_mixture"],
                        "delta_loo_vs_matched_independence": (
                            scores["loo_model_average"] - scores["matched_conditional_independence"]
                        ),
                    })

    expected = len(F0_RANGE_NAMES) * len(TRUTH_ANCHORS) * len(sample_sizes) * len(seed_roots)
    if len(rows) != expected:
        raise ArithmeticError(f"expected {expected} cohort rows, got {len(rows)}")

    group_summaries = []
    groups = sorted({(row["f0_range"], row["truth_label"], row["n"]) for row in rows})
    for f0_name, truth_label, n in groups:
        group = [
            row for row in rows
            if (row["f0_range"], row["truth_label"], row["n"]) == (f0_name, truth_label, n)
        ]
        methods = sorted(group[0]["scores"])
        scores = {
            method: {
                "mean": float(np.mean([row["scores"][method] for row in group])),
                "sd_across_cohort_replicates": float(np.std([row["scores"][method] for row in group], ddof=1))
                if len(group) > 1 else None,
            }
            for method in methods
        }
        group_summaries.append({
            "f0_range": f0_name,
            "truth_label": truth_label,
            "truth_parameters": group[0]["truth_parameters"],
            "primary_informative_truth": group[0]["primary_informative_truth"],
            "n": int(n),
            "n_cohort_replicates": len(group),
            "scores": scores,
            "mean_delta_loo_vs_uniform": float(np.mean([row["delta_loo_vs_uniform"] for row in group])),
            "mean_delta_loo_vs_matched_independence": float(
                np.mean([row["delta_loo_vs_matched_independence"] for row in group])
            ),
        })

    protocol_run = (
        seed_roots == ROOT_SEEDS
        and sample_sizes == SAMPLE_SIZES
        and int(n_f0_nodes) == N_F0_NODES
        and check_registered_seeds
    )
    return {
        "record_type": "synthetic_development_signal_prediction",
        "scope": "synthetic paired baseline-signal prediction only; no treatment outcomes or patient data",
        "protocol_run": bool(protocol_run),
        "protocol_sha256": pre_run_protocol_hash,
        "code_sha256": pre_run_code_hashes,
        "protocol": {
            "scenario_seed_roots": [int(seed_roots[0]), int(seed_roots[-1])],
            "root_seed_count": len(seed_roots),
            "f0_ranges": {name: list(map(float, F0_RANGES[name])) for name in F0_RANGE_NAMES},
            "truth_anchors": [list(map(float, row)) for row in TRUTH_ANCHORS],
            "primary_truth_anchors": [list(map(float, row)) for row in PRIMARY_TRUTHS],
            "sample_sizes": list(sample_sizes),
            "candidate_grid": [list(map(float, row)) for row in _RECOVERY.parameter_grid()],
            "candidate_count": len(MODEL_GRID),
            "quadrature_nodes": int(n_f0_nodes),
            "nominal_parameters": [0.8, 0.5, 0.4, 0.6],
            "seed_root_check": "stage2_dev" if check_registered_seeds else "disabled (test only)",
            "primary_score": "cohort mean per-patient natural-log predictive density",
            "primary_contrasts": ["loo_model_average - uniform_grid_mixture", "loo_model_average - matched_conditional_independence"],
            "uncertainty_unit": "root seed; each root contributes the equal mean over 3 informative truths x 3 f0 ranges",
        },
        "primary_gate": _primary_gate(rows, sample_sizes, seed_roots, protocol_run=protocol_run),
        "group_summaries": group_summaries,
        "all_cohort_scores": rows,
    }


def write_diagnostic(path: Path | None = None, **kwargs):
    path = ROOT / "output" / "trackc_loo_predictive_check" / "diagnostic.json" if path is None else Path(path)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing diagnostic: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    result = build_diagnostic(**kwargs)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path, result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true", help="run the frozen synthetic diagnostic")
    parser.add_argument("--output", type=Path, default=None, help="optional JSON output path")
    args = parser.parse_args(argv)
    if not args.run:
        print("No run: pass --run to execute the frozen synthetic diagnostic.")
        return 0
    path, result = write_diagnostic(args.output)
    print(json.dumps({
        "output": str(path),
        "protocol_run": result["protocol_run"],
        "cohorts": len(result["all_cohort_scores"]),
        "primary_gate": result["primary_gate"]["status"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
