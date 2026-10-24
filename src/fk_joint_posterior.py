"""Joint empirical-Bayes grid utilities for Track B log growth parameters.

These primitives preserve rho-D dependence and include quadrature weights. The caller must fit the prior
from outer-training patients only and supply the per-case likelihood; this module does not fit patient data.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import logsumexp


def _axis_quadrature_weights(axis) -> np.ndarray:
    axis = np.asarray(axis, dtype=float)
    if axis.ndim != 1 or len(axis) < 2 or not np.isfinite(axis).all() or np.any(np.diff(axis) <= 0):
        raise ValueError("quadrature axis must be finite, strictly increasing, and contain at least two points")
    gaps = np.diff(axis)
    weights = np.empty_like(axis)
    weights[0] = gaps[0] / 2.0
    weights[-1] = gaps[-1] / 2.0
    if len(axis) > 2:
        weights[1:-1] = (gaps[:-1] + gaps[1:]) / 2.0
    return weights


def _grid_from_axes(rho_log_grid, D_log_grid) -> dict:
    rho_log_grid = np.asarray(rho_log_grid, dtype=float)
    D_log_grid = np.asarray(D_log_grid, dtype=float)
    rho_weights = _axis_quadrature_weights(rho_log_grid)
    D_weights = _axis_quadrature_weights(D_log_grid)
    rho_mesh, D_mesh = np.meshgrid(rho_log_grid, D_log_grid, indexing="ij")
    points = np.column_stack((rho_mesh.ravel(), D_mesh.ravel()))
    quadrature_weights = np.outer(rho_weights, D_weights).ravel()
    return {"rho_log_grid": rho_log_grid, "D_log_grid": D_log_grid, "points": points,
            "quadrature_weights": quadrature_weights}


def make_log_parameter_grid(*, rho_bounds=(1e-5, 1e-1), D_bounds=(1e-5, 10.0),
                            n_coarse: int = 21) -> dict:
    """Create the coarse rho-D grid, equally spaced in log parameter coordinates."""
    if isinstance(n_coarse, (bool, np.bool_)) or int(n_coarse) != n_coarse or n_coarse < 2:
        raise ValueError("n_coarse must be an integer of at least two")
    axes = []
    for bounds, name in ((rho_bounds, "rho_bounds"), (D_bounds, "D_bounds")):
        values = np.asarray(bounds, dtype=float)
        if values.shape != (2,) or not np.isfinite(values).all() or values[0] <= 0 or values[0] >= values[1]:
            raise ValueError(f"{name} must contain two positive increasing finite values")
        axes.append(np.linspace(np.log(values[0]), np.log(values[1]), int(n_coarse)))
    grid = _grid_from_axes(*axes)
    grid["rho_bounds"] = tuple(float(value) for value in rho_bounds)
    grid["D_bounds"] = tuple(float(value) for value in D_bounds)
    grid["grid_type"] = "coarse log grid"
    return grid


def fit_log_parameter_prior(training_log_parameter_pairs, *, covariance_floor: float = 1e-6) -> dict:
    """Estimate a joint Gaussian prior from outer-training (log rho, log D) pairs."""
    pairs = np.asarray(training_log_parameter_pairs, dtype=float)
    covariance_floor = float(covariance_floor)
    if pairs.ndim != 2 or pairs.shape[1] != 2 or len(pairs) < 2 or not np.isfinite(pairs).all():
        raise ValueError("training_log_parameter_pairs must be finite with shape (n>=2, 2)")
    if not np.isfinite(covariance_floor) or covariance_floor <= 0:
        raise ValueError("covariance_floor must be positive and finite")
    mean = pairs.mean(axis=0)
    covariance = np.cov(pairs, rowvar=False, ddof=1) + covariance_floor * np.eye(2)
    if not np.isfinite(covariance).all() or np.linalg.eigvalsh(covariance).min() <= 0:
        raise ValueError("training parameters did not yield a positive-definite covariance")
    return {"mean": mean, "covariance": covariance, "n_training_pairs": int(len(pairs)),
            "parameter_order": ("log_rho", "log_D"), "covariance_floor": covariance_floor}


def log_gaussian_prior_density(points, prior: dict) -> np.ndarray:
    """Evaluate a bivariate Gaussian log density on joint log-parameter points."""
    points = np.asarray(points, dtype=float)
    mean = np.asarray(prior["mean"], dtype=float)
    covariance = np.asarray(prior["covariance"], dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError("points must be a finite array with shape (n, 2)")
    if mean.shape != (2,) or covariance.shape != (2, 2) or not np.isfinite(mean).all() or not np.isfinite(covariance).all():
        raise ValueError("prior must contain a finite length-2 mean and 2x2 covariance")
    try:
        factor = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError("prior covariance must be positive definite") from exc
    centered = (points - mean).T
    whitened = np.linalg.solve(factor, centered)
    quadratic = np.sum(whitened**2, axis=0)
    log_determinant = 2.0 * np.log(np.diag(factor)).sum()
    return -0.5 * (2.0 * math.log(2.0 * math.pi) + log_determinant + quadratic)


def normalize_grid_posterior(log_likelihood, log_prior, quadrature_weights) -> dict:
    """Normalize likelihood × prior on a grid using its nonuniform cell-quadrature weights."""
    log_likelihood = np.asarray(log_likelihood, dtype=float)
    log_prior = np.asarray(log_prior, dtype=float)
    quadrature_weights = np.asarray(quadrature_weights, dtype=float)
    if (log_likelihood.ndim != 1 or log_prior.shape != log_likelihood.shape
            or quadrature_weights.shape != log_likelihood.shape or len(log_likelihood) == 0):
        raise ValueError("likelihood, prior, and quadrature weights must be matching nonempty 1D arrays")
    if (np.isnan(log_likelihood).any() or np.isnan(log_prior).any() or not np.isfinite(quadrature_weights).all()
            or np.any(quadrature_weights <= 0)):
        raise ValueError("log densities may not contain NaN and quadrature weights must be positive finite")
    log_mass = log_likelihood + log_prior + np.log(quadrature_weights)
    log_evidence = float(logsumexp(log_mass))
    if not np.isfinite(log_evidence):
        raise ValueError("all grid posterior masses are zero or nonfinite")
    probabilities = np.exp(log_mass - log_evidence)
    return {"probabilities": probabilities, "log_evidence": log_evidence,
            "normalization": "likelihood × joint prior × log-grid quadrature weight"}


def refine_log_parameter_grid(coarse_grid: dict, coarse_log_posterior, *,
                              mass_fraction: float = 0.9, n_refine: int = 11) -> dict:
    """Add local log-grid resolution over the smallest high-posterior-mass region available on the coarse grid."""
    log_posterior = np.asarray(coarse_log_posterior, dtype=float)
    points = np.asarray(coarse_grid["points"], dtype=float)
    if log_posterior.ndim != 1 or log_posterior.shape != (len(points),) or np.isnan(log_posterior).any():
        raise ValueError("coarse_log_posterior must be a non-NaN vector aligned to coarse-grid points")
    if not 0 < mass_fraction < 1 or not np.isfinite(mass_fraction):
        raise ValueError("mass_fraction must be strictly between zero and one")
    if isinstance(n_refine, (bool, np.bool_)) or int(n_refine) != n_refine or n_refine < 3:
        raise ValueError("n_refine must be an integer of at least three")

    posterior = normalize_grid_posterior(
        log_posterior, np.zeros_like(log_posterior), coarse_grid["quadrature_weights"]
    )["probabilities"]
    order = np.argsort(posterior)[::-1]
    cumulative = np.cumsum(posterior[order])
    selected = order[:int(np.searchsorted(cumulative, mass_fraction, side="left")) + 1]
    rho_axis = np.asarray(coarse_grid["rho_log_grid"], dtype=float)
    D_axis = np.asarray(coarse_grid["D_log_grid"], dtype=float)
    lower = points[selected].min(axis=0)
    upper = points[selected].max(axis=0)
    for dimension, axis in enumerate((rho_axis, D_axis)):
        if upper[dimension] - lower[dimension] < np.finfo(float).eps:
            center = lower[dimension]
            step = float(np.median(np.diff(axis)))
            lower[dimension], upper[dimension] = center - step, center + step
        else:
            padding = float(np.median(np.diff(axis))) / 2.0
            lower[dimension] -= padding
            upper[dimension] += padding
    lower = np.maximum(lower, [rho_axis[0], D_axis[0]])
    upper = np.minimum(upper, [rho_axis[-1], D_axis[-1]])
    fine_rho = np.linspace(lower[0], upper[0], int(n_refine))
    fine_D = np.linspace(lower[1], upper[1], int(n_refine))
    refined = _grid_from_axes(np.unique(np.concatenate((rho_axis, fine_rho))),
                              np.unique(np.concatenate((D_axis, fine_D))))
    refined["rho_bounds"] = coarse_grid["rho_bounds"]
    refined["D_bounds"] = coarse_grid["D_bounds"]
    refined["grid_type"] = "coarse plus local high-mass refinement"
    refined["refinement_mass_fraction"] = float(mass_fraction)
    refined["refinement_region_log"] = {"rho": (float(lower[0]), float(upper[0])),
                                        "D": (float(lower[1]), float(upper[1]))}
    return refined


def sample_joint_log_parameters(points, posterior_probabilities, *, n_draws: int = 1, seed: int = 0) -> np.ndarray:
    """Sample joint (log rho, log D) grid nodes without breaking their posterior dependence."""
    points = np.asarray(points, dtype=float)
    probabilities = np.asarray(posterior_probabilities, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError("points must be finite with shape (n, 2)")
    if (probabilities.ndim != 1 or probabilities.shape != (len(points),) or not np.isfinite(probabilities).all()
            or np.any(probabilities < 0) or probabilities.sum() <= 0):
        raise ValueError("posterior_probabilities must be nonnegative finite weights aligned to points")
    if isinstance(n_draws, (bool, np.bool_)) or int(n_draws) != n_draws or n_draws < 1:
        raise ValueError("n_draws must be a positive integer")
    probabilities = probabilities / probabilities.sum()
    indices = np.random.default_rng(seed).choice(len(points), size=int(n_draws), p=probabilities)
    return points[indices]
