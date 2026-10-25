"""Synthetic feasibility gate for the Track B growth-zone count/location law.

This module uses only generated arrays.  Its spatial field is an analytic
anisotropic reach proxy, not a numerical PDE solution or patient MRI.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Any

import numpy as np
from scipy.ndimage import gaussian_filter, shift
from scipy.optimize import minimize
from scipy.special import expit, gammaln
from scipy.stats import nbinom


GRID_SHAPE = (8, 8, 8)
N_VOXELS = int(np.prod(GRID_SHAPE))


def log_elementary_symmetric(scores: np.ndarray, k: int) -> float:
    """Return log(e_k(exp(scores))) using O(n*k) dynamic programming."""
    eta = np.asarray(scores, dtype=np.float64).reshape(-1)
    if not np.all(np.isfinite(eta)):
        raise ValueError("scores must be finite")
    if isinstance(k, (bool, np.bool_)) or int(k) != k or not 0 <= int(k) <= eta.size:
        raise ValueError("k must be an integer between zero and the score count")
    k = int(k)
    if k == 0:
        return 0.0
    dp = np.full(k + 1, -np.inf, dtype=np.float64)
    dp[0] = 0.0
    seen = 0
    for score in eta:
        seen += 1
        upper = min(k, seen)
        dp[1:upper + 1] = np.logaddexp(dp[1:upper + 1], score + dp[:upper])
    return float(dp[k])


def sample_conditional_bernoulli(scores: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    """Draw an exact-size weighted subset with probability proportional to exp(sum score)."""
    eta = np.asarray(scores, dtype=np.float64).reshape(-1)
    if not np.all(np.isfinite(eta)):
        raise ValueError("scores must be finite")
    if isinstance(k, (bool, np.bool_)) or int(k) != k or not 0 <= int(k) <= eta.size:
        raise ValueError("k must be an integer between zero and the score count")
    k = int(k)
    n = eta.size
    table = np.full((n + 1, k + 1), -np.inf, dtype=np.float64)
    table[:, 0] = 0.0
    for i, score in enumerate(eta, start=1):
        upper = min(i, k)
        table[i, 1:upper + 1] = np.logaddexp(
            table[i - 1, 1:upper + 1], score + table[i - 1, :upper]
        )
    selected = np.zeros(n, dtype=bool)
    remaining = k
    for i in range(n, 0, -1):
        if remaining == 0:
            break
        log_p = eta[i - 1] + table[i - 1, remaining - 1] - table[i, remaining]
        probability = float(np.clip(np.exp(log_p), 0.0, 1.0))
        if rng.random() < probability:
            selected[i - 1] = True
            remaining -= 1
    if remaining != 0:
        raise RuntimeError("conditional Bernoulli sampler failed to allocate requested count")
    return selected


def conditional_location_log_probability(mask: np.ndarray, scores: np.ndarray) -> float:
    """Log P(mask | its count, scores) under the fixed-cardinality law."""
    y = np.asarray(mask, dtype=bool).reshape(-1)
    eta = np.asarray(scores, dtype=np.float64).reshape(-1)
    if y.size != eta.size:
        raise ValueError("mask and score lengths differ")
    k = int(y.sum())
    return float(np.dot(y.astype(np.float64), eta) - log_elementary_symmetric(eta, k))


def solve_mass_intercept(scores: np.ndarray, expected_count: float) -> np.ndarray:
    """Return independent Bernoulli probabilities whose sum equals expected_count."""
    eta = np.asarray(scores, dtype=np.float64).reshape(-1)
    if eta.size == 0 or not np.all(np.isfinite(eta)):
        raise ValueError("scores must be a nonempty finite vector")
    target = float(expected_count)
    if not np.isfinite(target) or not 0 <= target <= eta.size:
        raise ValueError("expected_count must lie within the candidate set")
    if target == 0:
        return np.zeros_like(eta)
    if target == eta.size:
        return np.ones_like(eta)
    lo, hi = -60.0 - float(np.max(eta)), 60.0 - float(np.min(eta))
    for _ in range(100):
        mid = (lo + hi) / 2.0
        total = float(expit(mid + eta).sum())
        if total < target:
            lo = mid
        else:
            hi = mid
    return expit(((lo + hi) / 2.0) + eta)


def independent_mask_log_probability(mask: np.ndarray, probabilities: np.ndarray) -> float:
    y = np.asarray(mask, dtype=bool).reshape(-1)
    p = np.clip(np.asarray(probabilities, dtype=np.float64).reshape(-1), 1e-12, 1.0 - 1e-12)
    if y.size != p.size:
        raise ValueError("mask and probability lengths differ")
    return float(np.sum(np.where(y, np.log(p), np.log1p(-p))))


def negative_binomial_logpmf(count: np.ndarray | int, mean: np.ndarray | float, dispersion: float) -> np.ndarray:
    """NB log PMF parameterized by mean and Var(K)=mean+dispersion*mean^2."""
    k = np.asarray(count, dtype=np.float64)
    mu = np.asarray(mean, dtype=np.float64)
    alpha = float(dispersion)
    if alpha <= 0 or not np.isfinite(alpha) or np.any(mu <= 0) or not np.all(np.isfinite(mu)):
        raise ValueError("mean and dispersion must be positive and finite")
    size = 1.0 / alpha
    return (
        gammaln(k + size) - gammaln(size) - gammaln(k + 1.0)
        + size * (np.log(size) - np.log(size + mu))
        + k * (np.log(mu) - np.log(size + mu))
    )


def fit_nb_regression(features: np.ndarray, counts: np.ndarray) -> dict[str, Any]:
    """Fit a low-dimensional NB mean model and shared dispersion by MLE."""
    x = np.asarray(features, dtype=np.float64)
    y = np.asarray(counts, dtype=np.int64).reshape(-1)
    if x.ndim != 2 or x.shape[0] != y.size or x.shape[0] < 5 or x.shape[1] < 1:
        raise ValueError("features and counts have incompatible or insufficient rows")
    if not np.all(np.isfinite(x)) or np.any(y < 0):
        raise ValueError("features must be finite and counts nonnegative")

    mean_count = max(float(y.mean()), 1e-3)
    initial_beta = np.zeros(x.shape[1], dtype=np.float64)
    initial_beta[0] = np.log(mean_count)
    initial = np.concatenate([initial_beta, [np.log(0.4)]])

    def objective(parameters: np.ndarray) -> float:
        beta = parameters[:-1]
        alpha = float(np.exp(parameters[-1]))
        mu = np.exp(np.clip(x @ beta, -8.0, 8.0))
        return -float(np.sum(negative_binomial_logpmf(y, mu, alpha)))

    bounds = [(-8.0, 8.0)] * x.shape[1] + [(-5.0, 2.5)]
    result = minimize(objective, initial, method="L-BFGS-B", bounds=bounds)
    if not np.all(np.isfinite(result.x)):
        raise RuntimeError("NB regression returned non-finite parameters")
    return {
        "coefficients": result.x[:-1].copy(),
        "dispersion": float(np.exp(result.x[-1])),
        "success": bool(result.success),
        "objective": float(result.fun),
    }


def _nb_pmf_on_support(mean: float, dispersion: float, max_count: int) -> np.ndarray:
    if mean <= 0 or dispersion <= 0:
        raise ValueError("mean and dispersion must be positive")
    size = 1.0 / dispersion
    probability = size / (size + mean)
    pmf = nbinom.pmf(np.arange(max_count + 1), size, probability)
    total = float(pmf.sum())
    if not np.isfinite(total) or total <= 0:
        raise RuntimeError("NB mass on count support is invalid")
    return pmf / total


def _predict_count_distribution(model: dict[str, Any], features: np.ndarray, max_count: int) -> tuple[float, float, np.ndarray]:
    mean = float(np.exp(np.clip(np.asarray(features) @ model["coefficients"], -8.0, 8.0)))
    pmf = _nb_pmf_on_support(mean, float(model["dispersion"]), max_count)
    truncated_mean = float(np.dot(np.arange(max_count + 1), pmf))
    return mean, truncated_mean, pmf


def _count_features(previous_count: int, interval_ratio: float) -> np.ndarray:
    return np.array([1.0, np.log1p(previous_count), np.log(interval_ratio)], dtype=np.float64)


def _true_count_mean(previous_count: int, interval_ratio: float) -> float:
    return float(np.exp(0.9 + 0.7 * np.log1p(previous_count) + 0.5 * np.log(interval_ratio)))


def _reach_score(axis: int) -> np.ndarray:
    coordinate = np.stack(np.meshgrid(*[np.linspace(-1.0, 1.0, n) for n in GRID_SHAPE], indexing="ij"), axis=-1)
    scales = np.full(3, 0.45, dtype=np.float64)
    scales[axis] = 1.8
    return -np.sqrt(np.sum(coordinate**2 / scales, axis=-1)).reshape(-1)


def _standardize(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=np.float64)
    scale = float(x.std())
    return (x - float(x.mean())) / (scale if scale > 1e-8 else 1.0)


def _make_spatial_case(
    rng: np.random.Generator,
    *,
    count_dispersion: float,
    cluster_strength: float,
    registration_shift: int,
    pde_misspecified: bool,
) -> dict[str, Any]:
    previous_count = int(rng.negative_binomial(8.0, 8.0 / 16.0))
    interval_ratio = float(np.exp(rng.uniform(-0.35, 0.35)))
    mean = _true_count_mean(previous_count, interval_ratio)
    size = 1.0 / count_dispersion
    count_probability = size / (size + mean)
    count = int(rng.negative_binomial(size, count_probability))
    count = min(count, N_VOXELS)

    true_axis = int(rng.integers(0, 3))
    model_axis = true_axis
    if pde_misspecified:
        model_axis = (true_axis + int(rng.integers(1, 3))) % 3
    true_pde = _standardize(_reach_score(true_axis))
    model_pde = _standardize(_reach_score(model_axis))
    coordinates = np.stack(np.meshgrid(*[np.linspace(-1.0, 1.0, n) for n in GRID_SHAPE], indexing="ij"), axis=-1)
    perpendicular = [axis for axis in range(3) if axis != true_axis]
    tract = _standardize(-np.sqrt(np.sum(coordinates[..., perpendicular] ** 2, axis=-1)).reshape(-1))

    white_noise = rng.normal(size=GRID_SHAPE)
    correlated_field = gaussian_filter(white_noise, sigma=1.0, mode="reflect").reshape(-1)
    correlated_field = _standardize(correlated_field)
    true_score = 0.75 * true_pde + 0.65 * tract + cluster_strength * correlated_field
    predicted_score = 0.75 * model_pde + 0.65 * tract
    if registration_shift:
        predicted_score = shift(
            predicted_score.reshape(GRID_SHAPE),
            shift=(registration_shift, 0, 0),
            order=1,
            mode="nearest",
        ).reshape(-1)
    mask = sample_conditional_bernoulli(true_score, count, rng)
    return {
        "previous_count": previous_count,
        "interval_ratio": interval_ratio,
        "count": count,
        "mask": mask,
        "scores": predicted_score,
    }


def _average_precision(mask: np.ndarray, scores: np.ndarray) -> float:
    y = np.asarray(mask, dtype=bool).reshape(-1)
    s = np.asarray(scores, dtype=np.float64).reshape(-1)
    positives = int(y.sum())
    if positives == 0:
        return float("nan")
    order = np.argsort(-s, kind="mergesort")
    ranked = y[order].astype(np.float64)
    precision_at_rank = np.cumsum(ranked) / (np.arange(ranked.size) + 1.0)
    return float(np.sum(precision_at_rank * ranked) / positives)


def _summarize_cell(rows: list[dict[str, float]]) -> dict[str, float]:
    keys = rows[0].keys()
    summary: dict[str, float] = {}
    for key in keys:
        values = np.asarray([row[key] for row in rows], dtype=np.float64)
        finite = np.isfinite(values)
        summary[key] = float(np.mean(values[finite])) if finite.any() else float("nan")
    return summary


def run_synthetic_gate(
    *,
    seeds: tuple[int, ...] = (20261025, 20261026, 20261027),
    n_train: int = 120,
    n_test: int = 50,
) -> dict[str, Any]:
    """Compare global, mean-count, and joint count/location laws on generated cases."""
    if n_train < 20 or n_test < 20:
        raise ValueError("use at least 20 synthetic training and test cases per condition")
    condition_rows: dict[tuple[float, float, int, int], list[dict[str, float]]] = {}
    count_cells = list(product((0.2, 1.0), (0.0, 1.2), (0, 1), (False, True)))

    for seed, cell in product(seeds, count_cells):
        dispersion, cluster, reg_shift, pde_wrong = cell
        seed_sequence = np.random.SeedSequence([int(seed), int(round(dispersion * 10)), int(round(cluster * 10)), reg_shift, int(pde_wrong)])
        rng = np.random.default_rng(seed_sequence)

        training_previous = rng.negative_binomial(8.0, 8.0 / 16.0, size=n_train)
        training_ratio = np.exp(rng.uniform(-0.35, 0.35, size=n_train))
        training_counts = np.array([
            rng.negative_binomial(
                1.0 / dispersion,
                (1.0 / dispersion) / ((1.0 / dispersion) + _true_count_mean(int(previous), float(ratio))),
            )
            for previous, ratio in zip(training_previous, training_ratio)
        ], dtype=np.int64)
        training_counts = np.minimum(training_counts, N_VOXELS)
        history_model = fit_nb_regression(
            np.vstack([_count_features(int(p), float(r)) for p, r in zip(training_previous, training_ratio)]),
            training_counts,
        )
        global_model = fit_nb_regression(np.ones((n_train, 1)), training_counts)

        key = (float(dispersion), float(cluster), int(reg_shift), int(pde_wrong))
        condition_rows.setdefault(key, [])
        for _ in range(n_test):
            case = _make_spatial_case(
                rng,
                count_dispersion=dispersion,
                cluster_strength=cluster,
                registration_shift=reg_shift,
                pde_misspecified=pde_wrong,
            )
            history_features = _count_features(case["previous_count"], case["interval_ratio"])
            _, history_mean, history_pmf = _predict_count_distribution(history_model, history_features, N_VOXELS)
            _, global_mean, global_pmf = _predict_count_distribution(global_model, np.ones(1), N_VOXELS)
            observed_count = int(case["count"])
            mask = case["mask"]
            scores = case["scores"]

            log_location = conditional_location_log_probability(mask, scores)
            count_nll_history = -float(np.log(max(history_pmf[observed_count], 1e-300)))
            count_nll_global = -float(np.log(max(global_pmf[observed_count], 1e-300)))
            joint_nll_history = count_nll_history - log_location
            joint_nll_global = count_nll_global - log_location

            p_history = solve_mass_intercept(scores, history_mean)
            p_global = solve_mass_intercept(scores, global_mean)
            mass_nll_history = -independent_mask_log_probability(mask, p_history)
            mass_nll_global = -independent_mask_log_probability(mask, p_global)

            # The NB intervals are truncated to the feasible count support.
            history_cdf = np.cumsum(history_pmf)
            global_cdf = np.cumsum(global_pmf)
            history_lo90 = int(np.searchsorted(history_cdf, 0.05, side="left"))
            history_hi90 = int(np.searchsorted(history_cdf, 0.95, side="left"))
            global_lo90 = int(np.searchsorted(global_cdf, 0.05, side="left"))
            global_hi90 = int(np.searchsorted(global_cdf, 0.95, side="left"))
            ap_raw = _average_precision(mask, scores)
            ap_mass = _average_precision(mask, np.log(np.clip(p_history, 1e-12, 1.0 - 1e-12) / np.clip(1.0 - p_history, 1e-12, 1.0)))

            condition_rows[key].append({
                "count_mae_history": abs(history_mean - observed_count),
                "count_mae_global": abs(global_mean - observed_count),
                "count_nll_history": count_nll_history,
                "count_nll_global": count_nll_global,
                "coverage90_history": float(history_lo90 <= observed_count <= history_hi90),
                "coverage90_global": float(global_lo90 <= observed_count <= global_hi90),
                "conditional_location_nll": -log_location,
                "joint_nll_history": joint_nll_history,
                "joint_nll_global": joint_nll_global,
                "joint_nll_expected_mass_history": mass_nll_history,
                "joint_nll_expected_mass_global": mass_nll_global,
                "history_joint_minus_expected_mass_nll": joint_nll_history - mass_nll_history,
                "history_joint_minus_global_joint_nll": joint_nll_history - joint_nll_global,
                "ap_raw": ap_raw,
                "ap_after_mass_intercept": ap_mass,
                "ap_invariance_abs_error": abs(ap_raw - ap_mass),
                "ap_evaluable": float(np.isfinite(ap_raw)),
                "location_nll_gain_vs_uniform": math_log_choose(N_VOXELS, observed_count) + log_location,
            })

    summarized = []
    for key, rows in sorted(condition_rows.items()):
        dispersion, cluster, reg_shift, pde_wrong = key
        row = _summarize_cell(rows)
        row.update({
            "count_dispersion": dispersion,
            "spatial_cluster_strength": cluster,
            "registration_shift_voxels": reg_shift,
            "pde_misspecified": bool(pde_wrong),
            "n_test": len(rows),
            "n_ap_evaluable": int(round(row["ap_evaluable"] * len(rows))),
        })
        summarized.append(row)

    return {
        "scope": "synthetic-only; no patient data loaded or scored",
        "spatial_generator": "8x8x8 analytic anisotropic reach score plus optional smooth latent field; not a numerical PDE",
        "count_generator": "negative-binomial next-count model conditional on previous count and interval ratio",
        "count_location_law": "truncated NB count times conditional-Bernoulli fixed-count location law",
        "seeds": list(seeds),
        "n_train_per_seed_condition": n_train,
        "n_test_per_seed_condition": n_test,
        "n_conditions": len(count_cells),
        "conditions": summarized,
        "interpretation_limits": [
            "This is a controlled implementation and falsification exercise, not evidence about MRI performance.",
            "The true spatial process is generated from the same broad score family plus an optional smooth latent field; it is not biological ground truth.",
            "The per-voxel probabilities after a scalar mass intercept preserve the score ranking by construction.",
            "Conditional Bernoulli adds a fixed-count law but no explicit attractive neighbor interaction.",
        ],
    }


def math_log_choose(n: int, k: int) -> float:
    return float(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("output/track_b_growth_zone_joint_law_synthetic.json"))
    parser.add_argument("--train", type=int, default=120)
    parser.add_argument("--test", type=int, default=50)
    parser.add_argument("--seeds", type=int, nargs="+", default=[20261025, 20261026, 20261027])
    args = parser.parse_args()
    result = run_synthetic_gate(seeds=tuple(args.seeds), n_train=args.train, n_test=args.test)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "conditions": len(result["conditions"]),
        "total_test_cases": len(result["conditions"]) * len(args.seeds) * args.test,
        "scope": result["scope"],
    }, indent=2))
