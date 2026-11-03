import numpy as np
import pytest

from src import fk_two_history_shape as th
from src.fk_hybrid_forecast import signed_distance_volume_baseline_draws

COV = {"sigma_b_mm": 1.0, "tau_b_mm": 0.3, "ell_mm": 4.0}
SHAPE = (24, 24, 24)


def ball(c, r, shape=SHAPE):
    g = np.indices(shape)
    return sum((g[i] - c[i]) ** 2 for i in range(3)) <= r * r


@pytest.fixture
def setup():
    h1, h2 = ball((11, 12, 12), 4), ball((12, 12, 12), 5)
    dom = np.ones(SHAPE, bool)
    dom[:2] = dom[-2:] = False
    return h1, h2, dom, np.array([300, 420, 520], dtype=np.int64)


def run(setup, **kw):
    h1, h2, dom, k = setup
    return th.two_history_shape_draws(h1, h2, 30, 30, k, domain=dom, boundary_covariance=COV,
                                      n_anchors=64, **kw)


def test_exact_volume_and_support(setup):
    _, _, dom, k = setup
    r = run(setup)
    assert np.array_equal(r["mask_draws"].reshape(3, -1).sum(1), k)
    assert np.array_equal(r["pre_boundary_mask_draws"].reshape(3, -1).sum(1), k)
    assert not (r["mask_draws"] & ~dom).any() and not (r["pre_boundary_mask_draws"] & ~dom).any()


def test_seed_determinism_and_shared_seeds(setup):
    a, b = run(setup, seed=3), run(setup, seed=3)
    assert np.array_equal(a["mask_draws"], b["mask_draws"])
    c = run(setup, boundary_seeds=th.candidate_boundary_seeds(3, 3))
    assert np.array_equal(a["mask_draws"], c["mask_draws"])
    d = run(setup, boundary_seeds=th.candidate_boundary_seeds(4, 3))
    assert not np.array_equal(a["mask_draws"], d["mask_draws"])
    with pytest.raises(ValueError):
        run(setup, boundary_seeds=[1, 2])


def test_candidate_seed_derivation_matches_hybrid_source():
    n = 5
    _, bs = np.random.SeedSequence(9).spawn(2)
    expect = [int(s.generate_state(1)[0]) for s in bs.spawn(n)]
    assert th.candidate_boundary_seeds(9, n) == expect


def test_reduces_to_last_mask_rule_when_h1_equals_h2(setup):
    _, h2, dom, k = setup
    r = th.two_history_shape_draws(h2, h2, 30, 30, k, domain=dom, apply_boundary=False)
    last = signed_distance_volume_baseline_draws(
        h2, k, domain=dom, boundary_covariance={"sigma_b_mm": 0.0, "tau_b_mm": 0.0, "ell_mm": 4.0})
    assert np.array_equal(r["pre_boundary_mask_draws"], last)


def test_differs_from_last_mask_when_history_moves(setup):
    h1, h2, dom, k = setup
    moved = th.two_history_shape_draws(h1, h2, 30, 60, k, domain=dom, apply_boundary=False)
    still = th.two_history_shape_draws(h2, h2, 30, 60, k, domain=dom, apply_boundary=False)
    assert not np.array_equal(moved["pre_boundary_mask_draws"], still["pre_boundary_mask_draws"])
    assert moved["factor"] > 0


def test_tiny_and_zero_gap_fall_back(setup):
    h1, h2, dom, k = setup
    still = th.two_history_shape_draws(h2, h2, 30, 30, k, domain=dom, apply_boundary=False)
    for gap in (0.0, 0.5):
        r = th.two_history_shape_draws(h1, h2, gap, 30, k, domain=dom, apply_boundary=False)
        assert r["factor"] == 0.0 and r["tiny_gap_fallback"]
        assert np.array_equal(r["pre_boundary_mask_draws"], still["pre_boundary_mask_draws"])


def test_factor_cap_and_shrinkage():
    assert th.extrapolation_factor(1.0, 1e6) == th.MAX_FACTOR
    assert th.extrapolation_factor(7.0, 7.0) == pytest.approx(0.5)  # w = 7/14
    assert th.extrapolation_factor(0.2, 10) == 0.0
    with pytest.raises(ValueError):
        th.extrapolation_factor(10, 0)


def test_empty_masks(setup):
    h1, h2, dom, k = setup
    z = np.zeros(SHAPE, bool)
    r = th.two_history_shape_draws(z, h2, 30, 30, k, domain=dom, apply_boundary=False)
    assert r["h1_unusable_fallback"] and r["factor"] == 0.0
    with pytest.raises(ValueError):
        th.two_history_shape_draws(h1, z, 30, 30, k, domain=dom, apply_boundary=False)


def test_shrinking_history_still_matches_volume(setup):
    h1, h2, dom, k = setup
    r = th.two_history_shape_draws(h2, h1, 30, 30, k, domain=dom, boundary_covariance=COV, n_anchors=64)
    assert np.array_equal(r["mask_draws"].reshape(3, -1).sum(1), k)
    assert r["factor"] > 0


def test_pre_boundary_path_deterministic_and_needs_no_covariance(setup):
    h1, h2, dom, k = setup
    a = th.two_history_shape_draws(h1, h2, 30, 30, k, domain=dom, apply_boundary=False)
    b = th.two_history_shape_draws(h1, h2, 30, 30, k, domain=dom, apply_boundary=False, seed=99)
    assert np.array_equal(a["mask_draws"], b["mask_draws"]) and not a["boundary_applied"]
    with pytest.raises(ValueError):
        th.two_history_shape_draws(h1, h2, 30, 30, k, domain=dom)


def test_input_validation(setup):
    h1, h2, dom, k = setup
    with pytest.raises(ValueError):
        th.two_history_shape_draws(h1, h2, 30, 30, np.array([dom.sum()]), domain=dom, apply_boundary=False)
    with pytest.raises(ValueError):
        th.two_history_shape_draws(h1, h2, 30, 30, k, domain=h2 & False, apply_boundary=False)


def test_consensus_mask(setup):
    h1, h2, dom, k = setup
    r = th.two_history_shape_draws(h1, h2, 30, 30, k, domain=dom, apply_boundary=False)
    m, p = th.consensus_mask(r["mask_draws"], dom)
    assert m.shape == SHAPE and p.max() <= 1.0 and not (m & ~dom).any()
