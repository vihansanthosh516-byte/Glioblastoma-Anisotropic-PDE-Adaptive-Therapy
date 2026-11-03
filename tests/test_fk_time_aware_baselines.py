"""Synthetic tests for the PROPOSED, NOT ADOPTED time-aware baselines. No patient data."""
import numpy as np
import pytest
import torch

from src import fk_time_aware_baselines as ta
from src.fk_learned_baselines import UNet3D

WIDTHS = (4, 8)  # small backbone keeps the tests fast; the structure is the same as the frozen one


def _pair():
    torch.manual_seed(20261102)
    frozen = UNet3D(WIDTHS)
    film = ta.FiLMUNet3D(WIDTHS)
    film.backbone.load_state_dict(frozen.state_dict(), strict=True)
    return frozen, film


def _bytes(t):
    return t.detach().contiguous().numpy().tobytes()


def test_film_layout_matches_frozen_unet():
    m = ta.FiLMUNet3D()
    assert m.norm_channels == [16, 16, 32, 32, 64, 64, 128, 128, 64, 64, 32, 32, 16, 16]
    assert m.film[-1].out_features == 1408
    assert sum(p.numel() for p in m.film.parameters()) == 23984


def test_zero_init_equals_frozen_outputs_and_backbone_gradients():
    frozen, film = _pair()
    rng = np.random.default_rng(20261105)
    x = torch.from_numpy(rng.random((1, 3, 16, 16, 16), dtype=np.float32))
    y = torch.from_numpy((rng.random((1, 1, 16, 16, 16)) > 0.5).astype(np.float32))
    from src.fk_learned_baselines import _loss
    prior = torch.are_deterministic_algorithms_enabled()
    torch.use_deterministic_algorithms(True)
    try:
        frozen.zero_grad()
        lf = _loss(frozen(x), y, 1.0)
        lf.backward()
        ref_out = frozen(x).detach()
        ref_grads = {n: p.grad.clone() for n, p in frozen.named_parameters()}
        for z in ((13, 13), (15, 155), (11, 1)):
            zt = ta.interval_features(*z)[None]
            film.zero_grad()
            out = film(x, zt)
            assert torch.equal(out, ref_out) and _bytes(out) == _bytes(ref_out)
            lz = _loss(out, y, 1.0)
            assert torch.equal(lz, lf.detach()) and _bytes(lz) == _bytes(lf)
            lz.backward()
            for n, p in film.backbone.named_parameters():
                assert torch.equal(p.grad, ref_grads[n]) and _bytes(p.grad) == _bytes(ref_grads[n]), n
            assert film.film[-1].weight.grad.abs().sum() > 0
    finally:
        torch.use_deterministic_algorithms(prior)


def test_interval_changes_output_with_fixed_masks():
    _, film = _pair()
    with torch.no_grad():
        torch.manual_seed(7)
        film.film[-1].weight.normal_(0, 0.1)
        film.film[-1].bias.zero_()
    x = torch.from_numpy(np.random.default_rng(1).random((1, 3, 16, 16, 16), dtype=np.float32))
    with torch.no_grad():
        base = film(x, ta.interval_features(13, 13)[None])
        fwd = film(x, ta.interval_features(15, 13)[None])
        hist = film(x, ta.interval_features(13, 155)[None])
        again = film(x, ta.interval_features(13, 13)[None])
    assert torch.equal(base, again)
    assert (fwd - base).abs().max() > 1e-4 and (hist - base).abs().max() > 1e-4


def test_rate_grid_and_cap():
    assert ta.RATE_GRID_MM_PER_DAY[0] == 0.0 and ta.RATE_GRID_MM_PER_DAY[-1] == 0.19
    assert len(ta.RATE_GRID_MM_PER_DAY) == 20
    assert ta.timed_radius_mm(0.19, 15) == pytest.approx(19.95) and ta.timed_radius_mm(0.19, 15) <= 20.0
    with pytest.raises(ValueError, match="cap"):
        ta.select_rate([], rate_grid=(0.0, 0.20))
    with pytest.raises(ValueError):
        ta.timed_radius_mm(0.1, 16)


def _ball(r_mm, shape=(40, 40, 40), h=2.0):
    g = np.ogrid[tuple(slice(0, s) for s in shape)]
    return sum(((g[i] - 20) * h) ** 2 for i in range(3)) <= r_mm ** 2


def _case(pid, gap, k_true, r0=8.0):
    target_r = r0 + k_true * 7 * gap
    return ta.TimedGrowthCase(pid, _ball(r0), np.ones((40, 40, 40), bool), (2.0, 2.0, 2.0), gap, _ball(target_r))


def test_selection_recovers_rate_and_ties_go_smaller():
    sel = ta.select_rate([_case("a", 11, 0.08), _case("b", 15, 0.08)])
    assert sel["rate_mm_per_day"] == pytest.approx(0.08, abs=0.011)
    # identical targets for every k (target = last core inside a 1-voxel support) -> all tie -> smallest k
    tiny = np.zeros((40, 40, 40), bool)
    tiny[20, 20, 20] = True
    tie = ta.TimedGrowthCase("t", tiny, tiny, (2.0, 2.0, 2.0), 13, tiny)
    sel = ta.select_rate([tie])
    assert sel["rate_mm_per_day"] == 0.0 and sel["equals_copy_last"]


def test_outer_selection_uses_training_patients_only():
    cases = [_case("p0", 12, 0.04), _case("p1", 13, 0.04), _case("p2", 14, 0.15), _case("p3", 11, 0.15)]
    folds = {"p0": 0, "p1": 0, "p2": 1, "p3": 1}
    out = ta.run_outer_timed_growth(cases, folds)
    assert out[0]["training_ids"] == ["p2", "p3"] and out[1]["training_ids"] == ["p0", "p1"]
    # fold 0 is fitted on fast growers only, fold 1 on slow growers only
    assert out[0]["selection"]["rate_mm_per_day"] > out[1]["selection"]["rate_mm_per_day"]
    assert set(out[0]["forecasts"]) == {"p0", "p1"}
