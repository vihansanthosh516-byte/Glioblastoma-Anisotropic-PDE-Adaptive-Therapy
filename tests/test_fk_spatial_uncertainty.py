from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.fk_spatial_uncertainty import (  # noqa: E402
    correlated_boundary_perturbation,
    deterministic_surface_anchors,
    energy_score_mask,
    exponential_covariance_matrix,
    fixed_boundary_score_pairs,
    marginal_brier_score,
    sample_anchor_displacements,
    spatial_variogram_score,
    volume_crps,
)


def ball(shape=(15, 15, 15), radius=4.0):
    coords = np.ogrid[tuple(slice(0, n) for n in shape)]
    center = tuple((n - 1) / 2 for n in shape)
    return sum((axis - center_i) ** 2 for axis, center_i in zip(coords, center)) <= radius**2


def test_surface_anchors_are_deterministic_and_on_the_mask_boundary():
    mask = ball()
    first = deterministic_surface_anchors(mask, spacing_mm=1.0, n_anchors=24)
    second = deterministic_surface_anchors(mask, spacing_mm=1.0, n_anchors=24)

    assert first.shape == (24, 3)
    assert np.array_equal(first, second)
    assert len(np.unique(first, axis=0)) == 24
    voxels = np.rint(first).astype(int)
    assert np.all(mask[tuple(voxels.T)])


def test_exponential_covariance_matches_the_declared_kernel_and_seed_is_reproducible():
    points = np.asarray([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    covariance = exponential_covariance_matrix(points, sigma_b_mm=2.0, tau_b_mm=0.5, ell_mm=5.0)

    assert np.isclose(covariance[0, 0], 4.25)
    assert np.isclose(covariance[0, 1], 4.0 * np.exp(-0.4))
    assert np.allclose(covariance, covariance.T)
    first = sample_anchor_displacements(points, sigma_b_mm=2.0, tau_b_mm=0.5, ell_mm=5.0, seed=12)
    second = sample_anchor_displacements(points, sigma_b_mm=2.0, tau_b_mm=0.5, ell_mm=5.0, seed=12)
    assert np.array_equal(first, second)


def test_correlated_boundary_perturbation_is_reproducible_and_zero_noise_is_identity():
    mask = ball()
    points = deterministic_surface_anchors(mask, n_anchors=32)
    assert np.array_equal(
        correlated_boundary_perturbation(mask, sigma_b_mm=0.0, tau_b_mm=0.0, ell_mm=5.0,
                                         spacing_mm=1.0, seed=3),
        mask,
    )
    first = correlated_boundary_perturbation(
        mask, sigma_b_mm=1.0, tau_b_mm=0.2, ell_mm=5.0, spacing_mm=1.0, seed=7, anchor_points_mm=points
    )
    second = correlated_boundary_perturbation(
        mask, sigma_b_mm=1.0, tau_b_mm=0.2, ell_mm=5.0, spacing_mm=1.0, seed=7, anchor_points_mm=points
    )
    assert np.array_equal(first, second)
    assert first.shape == mask.shape and first.any() and not first.all()


def test_fixed_boundary_pair_design_is_deterministic_and_within_distance_bins():
    mask = ball()
    first = fixed_boundary_score_pairs(mask, spacing_mm=1.0, n_voxels=64, band_mm=5.0)
    second = fixed_boundary_score_pairs(mask, spacing_mm=1.0, n_voxels=64, band_mm=5.0)

    assert np.array_equal(first["voxel_indices"], second["voxel_indices"])
    assert np.array_equal(first["pairs"], second["pairs"])
    assert np.all(first["distance_mm"] > 0)
    assert np.all(first["distance_mm"] <= 10.0)
    assert set(first["strata"]).issubset({0, 1, 2})


def test_spatial_score_design_is_anchored_to_history_before_target_is_seen():
    history = ball()
    draws = np.stack([history, np.roll(history, 1, axis=0)])
    target_a = np.roll(history, 1, axis=1)
    target_b = np.roll(history, 2, axis=2)
    design_a = fixed_boundary_score_pairs(history_mask=history, spacing_mm=1.0,
                                          n_voxels=64, band_mm=5.0)
    design_b = fixed_boundary_score_pairs(history_mask=history, spacing_mm=1.0,
                                          n_voxels=64, band_mm=5.0)

    assert np.array_equal(design_a["voxel_indices"], design_b["voxel_indices"])
    assert np.array_equal(design_a["pairs"], design_b["pairs"])
    assert spatial_variogram_score(draws, target_a, design_a)["score"] != \
        spatial_variogram_score(draws, target_b, design_b)["score"]


def test_variogram_brier_energy_and_volume_crps_match_hand_calculations():
    draws = np.asarray([
        [[[0, 0], [0, 1]]],
        [[[1, 0], [0, 0]]],
    ], dtype=bool)
    target = np.asarray([[[0, 1], [0, 0]]], dtype=bool)
    design = {
        "voxel_indices": np.asarray([0, 1, 2, 3]),
        "pairs": np.asarray([[0, 1], [2, 3]]),
        "distance_mm": np.asarray([1.0, 6.0]),
        "strata": np.asarray([0, 2]),
        "stratum_names": ("0-2 mm", "2-5 mm", "5-10 mm"),
    }

    variogram = spatial_variogram_score(draws, target, design)
    brier = marginal_brier_score(draws, target, design["voxel_indices"])
    energy = energy_score_mask(draws, target, design["voxel_indices"])

    assert np.isclose(variogram["score"], 0.25)
    assert variogram["strata"]["0-2 mm"]["n_pairs"] == 1
    assert variogram["strata"]["5-10 mm"]["n_pairs"] == 1
    assert np.isclose(brier, 0.375)
    assert np.isclose(energy, np.sqrt(2) - np.sqrt(2) / 4)
    assert np.isclose(volume_crps([0.0, 2.0], 1.0), 0.5)
