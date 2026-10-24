"""Synthetic Track B checks for anisotropy, boundary noise, and held-out forecasts."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src import fk_mask_inverse as fi  # noqa: E402
from src import fk_synthetic_validation as fv  # noqa: E402


def ball(shape=(24, 24, 24), radius=4.0):
    coords = np.ogrid[tuple(slice(0, n) for n in shape)]
    center = tuple((n - 1) / 2 for n in shape)
    return sum((axis - c) ** 2 for axis, c in zip(coords, center)) <= radius**2


def test_tensor_forward_matches_scalar_isotropic_forward():
    mask = ball(radius=3.0)
    domain = np.ones(mask.shape, bool)
    d6 = np.zeros(mask.shape + (6,), np.float32)
    d6[..., :3] = 0.08

    scalar = fi.forward_masks(mask, domain, 0.04, 0.08, [0.0, 8.0], 1.0)
    tensor = fi.forward_masks_tensor(mask, domain, 0.04, d6, [0.0, 8.0], 1.0)

    assert len(tensor) == 1
    assert np.array_equal(tensor[0], scalar[0])


def test_anisotropic_reference_expands_faster_along_high_diffusion_axis():
    mask = np.zeros((24, 24, 24), bool)
    mask[9:15, 9:15, 9:15] = True
    domain = np.ones(mask.shape, bool)
    d6 = np.zeros(mask.shape + (6,), np.float32)
    d6[..., 0], d6[..., 1], d6[..., 2] = 0.3, 0.03, 0.03

    pred = fi.forward_masks_tensor(mask, domain, 0.08, d6, [0.0, 20.0], 1.0)[0]
    extent = np.ptp(np.argwhere(pred), axis=0)

    assert extent[0] > extent[1]
    assert extent[0] > extent[2]


def test_signed_distance_boundary_noise_is_seeded_and_zero_is_identity():
    mask = ball(shape=(32, 32, 32), radius=7.0)
    zero = fv.perturb_boundary(mask, sigma_mm=0.0, correlation_mm=2.0, spacing_mm=1.0, seed=7)
    first = fv.perturb_boundary(mask, sigma_mm=0.8, correlation_mm=2.0, spacing_mm=1.0, seed=7)
    again = fv.perturb_boundary(mask, sigma_mm=0.8, correlation_mm=2.0, spacing_mm=1.0, seed=7)

    assert np.array_equal(zero, mask)
    assert first.dtype == bool and first.shape == mask.shape
    assert np.array_equal(first, again)
    assert not np.array_equal(first, mask)


def test_signed_distance_boundary_noise_rejects_unbounded_masks():
    with np.testing.assert_raises(ValueError):
        fv.perturb_boundary(np.zeros((8, 8, 8), bool), sigma_mm=1.0, correlation_mm=2.0, seed=1)


def test_volume_projection_hits_target_inside_domain():
    domain = ball(shape=(20, 20, 20), radius=8.0)
    mask = ball(shape=(20, 20, 20), radius=3.0) & domain
    target = int(mask.sum()) + 31

    grown = fv.project_mask_volume(mask, target, domain=domain, spacing_mm=2.0)
    shrunk = fv.project_mask_volume(grown, 50, domain=domain, spacing_mm=2.0)

    assert int(grown.sum()) == target
    assert int(shrunk.sum()) == 50
    assert np.all(grown <= domain) and np.all(shrunk <= domain)


def test_heldout_truth_does_not_change_fitted_forecast(monkeypatch):
    shape = (16, 16, 16)
    train = []
    for radius in (2.0, 2.5, 3.0):
        train.append(ball(shape=shape, radius=radius))
    days = [0.0, 10.0, 20.0]
    observed_fit_inputs = []

    def fake_fit(masks, fit_days, h, domain=None):
        observed_fit_inputs.append((len(masks), tuple(fit_days), [int(m.sum()) for m in masks]))
        return {"rho": 0.04, "D": 0.05, "sse_log_volume": 0.0, "dice": [1.0, 1.0]}

    monkeypatch.setattr(fi, "fit_fk_masks", fake_fit)
    truth_a = ball(shape=shape, radius=4.0)
    truth_b = np.roll(truth_a, 2, axis=0)

    result_a = fv.compare_next_scan(train, days, 30.0, truth_a, h=2.0)
    result_b = fv.compare_next_scan(train, days, 30.0, truth_b, h=2.0)

    assert observed_fit_inputs[0] == observed_fit_inputs[1]
    assert observed_fit_inputs[0][0:2] == (3, tuple(days))
    for arm in ("pde", "persistence", "volume_trend", "pde_shape_volume_trend"):
        assert np.array_equal(result_a["predictions"][arm], result_b["predictions"][arm])
        assert set(result_a["metrics"][arm]) >= {"dice", "hd95_mm", "spe_pct"}
    assert result_a["metrics"]["pde"] != result_b["metrics"]["pde"]


def test_next_scan_comparison_integrates_with_the_parameter_fitter():
    shape = (16, 16, 16)
    domain = np.ones(shape, bool)
    initial = ball(shape=shape, radius=2.5)
    days = [0.0, 20.0, 40.0, 60.0]
    reference = fi.forward_masks(initial, domain, 0.04, 0.05, days, 2.0)
    training_masks = [initial, *reference[:2]]

    result = fv.compare_next_scan(training_masks, days[:3], days[3], reference[2], h=2.0, domain=domain)

    assert result["fit"]["rho"] > 0 and result["fit"]["D"] > 0
    assert set(result["predictions"]) == {"pde", "persistence", "volume_trend", "pde_shape_volume_trend"}
    assert int(result["predictions"]["pde_shape_volume_trend"].sum()) == result["target_volume_trend_voxels"]
    assert all(mask.shape == shape for mask in result["predictions"].values())
    assert all(set(metrics) >= {"dice", "hd95_mm", "spe_pct"} for metrics in result["metrics"].values())
