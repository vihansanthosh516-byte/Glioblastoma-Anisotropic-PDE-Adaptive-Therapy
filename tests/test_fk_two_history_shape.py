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


def test_explicit_settings_preserve_defaults_and_control_factor(setup):
    default = th.TwoHistoryShapeSettings()
    assert default.max_factor == th.MAX_FACTOR
    assert default.min_history_gap_days == th.MIN_HISTORY_GAP_DAYS
    assert default.shrink_tau_days == th.SHRINK_TAU_DAYS
    assert default.max_sd_delta_mm == th.MAX_SD_DELTA_MM

    h1, h2, dom, k = setup
    custom = th.TwoHistoryShapeSettings(
        max_factor=0.25,
        min_history_gap_days=2.0,
        shrink_tau_days=0.0,
        max_sd_delta_mm=4.0,
    )
    assert th.extrapolation_factor(30.0, 60.0, settings=custom) == 0.25
    assert th.extrapolation_factor(1.0, 60.0, settings=custom) == 0.0
    tau_settings = th.TwoHistoryShapeSettings(
        max_factor=10.0,
        min_history_gap_days=0.0,
        shrink_tau_days=30.0,
    )
    assert th.extrapolation_factor(30.0, 30.0, settings=tau_settings) == 0.5

    result = th.two_history_shape_draws(
        h1, h2, 30.0, 60.0, k, domain=dom, apply_boundary=False, settings=custom
    )
    assert result["factor"] == 0.25
    assert result["settings"] == {
        "max_factor": 0.25,
        "min_history_gap_days": 2.0,
        "shrink_tau_days": 0.0,
        "max_sd_delta_mm": 4.0,
    }


def test_custom_distance_clip_changes_the_forecast():
    h1, h2 = ball((8, 12, 12), 4), ball((12, 12, 12), 5)
    low_clip = th.TwoHistoryShapeSettings(
        min_history_gap_days=0.0, shrink_tau_days=0.0, max_sd_delta_mm=1.0
    )
    high_clip = th.TwoHistoryShapeSettings(
        min_history_gap_days=0.0, shrink_tau_days=0.0, max_sd_delta_mm=20.0
    )
    low, _ = th.forecast_signed_distance(h1, h2, 30.0, 30.0, settings=low_clip)
    high, _ = th.forecast_signed_distance(h1, h2, 30.0, 30.0, settings=high_clip)
    assert not np.array_equal(low, high)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_factor": 0.0},
        {"max_factor": float("inf")},
        {"min_history_gap_days": -1.0},
        {"shrink_tau_days": -1.0},
        {"max_sd_delta_mm": 0.0},
        {"max_sd_delta_mm": float("nan")},
    ],
)
def test_explicit_settings_reject_invalid_values(kwargs):
    with pytest.raises(ValueError):
        th.TwoHistoryShapeSettings(**kwargs)


def test_zero_history_gap_is_last_mask_even_with_zero_threshold_and_tau():
    settings = th.TwoHistoryShapeSettings(min_history_gap_days=0.0, shrink_tau_days=0.0)
    assert th.extrapolation_factor(0.0, 30.0, settings=settings) == 0.0
    with pytest.raises(ValueError):
        th.extrapolation_factor(-1.0, 30.0, settings=settings)


def _inner_cv_fixture(n=20):
    outer_ids = [f"p{i:02d}" for i in range(n)]
    folds = []
    for fold in range(5):
        validation = outer_ids[fold * 4:(fold + 1) * 4]
        training = [pid for pid in outer_ids if pid not in validation]
        folds.append({"train_ids": training, "validation_ids": validation})
    return outer_ids, folds


TEST_IDS = ("outer-test",)


def test_inner_fold_plan_and_selector_are_training_only_with_frozen_ties():
    outer_ids, folds = _inner_cv_fixture()
    grid = {
        "f0": None,
        "loose": th.TwoHistoryShapeSettings(max_factor=4.0, shrink_tau_days=7.0),
        "damped": th.TwoHistoryShapeSettings(max_factor=1.0, shrink_tau_days=28.0),
    }
    for fold in folds:
        fold["scores"] = {
            "f0": {pid: 0.72 for pid in fold["validation_ids"]},
            "loose": {pid: 0.753 for pid in fold["validation_ids"]},
            "damped": {pid: 0.751 for pid in fold["validation_ids"]},
        }

    plan = th.validate_inner_fold_plan(outer_ids, folds, outer_test_ids=TEST_IDS)
    selected = th.select_two_history_settings(outer_ids, folds, grid, outer_test_ids=TEST_IDS)
    assert plan["n_outer_training_patients"] == 20
    assert plan["minimum_inner_training_patients"] == 10
    assert selected["selected_cell"] == "damped"
    assert selected["mean_inner_dice_by_cell"]["damped"] == pytest.approx(0.751)
    assert selected["inner_fold_selected_cells"] == ["damped"] * 5


def test_inner_fold_plan_rejects_test_ids_and_small_training_splits():
    outer_ids, folds = _inner_cv_fixture()
    folds[0]["validation_ids"].append("outer-test")
    with pytest.raises(ValueError, match="outside the outer-training set"):
        th.validate_inner_fold_plan(outer_ids, folds, outer_test_ids=TEST_IDS)

    outer_ids, folds = _inner_cv_fixture()
    with pytest.raises(ValueError, match="minimum inner-training size"):
        th.validate_inner_fold_plan(outer_ids, folds, outer_test_ids=TEST_IDS, minimum_inner_training_patients=17)

    outer_ids, folds = _inner_cv_fixture()
    grid = {"f0": None, "damped": th.TwoHistoryShapeSettings(max_factor=1.0, shrink_tau_days=28.0)}
    for fold in folds:
        fold["scores"] = {
            "f0": {pid: 0.7 for pid in fold["validation_ids"]},
            "damped": {pid: 0.7 for pid in fold["validation_ids"]},
        }
    folds[0]["scores"]["f0"]["outer-test"] = 1.0
    with pytest.raises(ValueError, match="score validation patients only"):
        th.select_two_history_settings(outer_ids, folds, grid, outer_test_ids=TEST_IDS)

    for fold in folds:
        fold["scores"]["f0"].pop("outer-test", None)
    selected = th.select_two_history_settings(outer_ids, folds, grid, outer_test_ids=TEST_IDS)
    assert selected["selected_cell"] == "f0"


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


def test_selector_rejects_outer_test_ids_in_training_set_and_requires_them():
    outer_ids, folds = _inner_cv_fixture()
    with pytest.raises(ValueError, match="overlap the outer-training set"):
        th.validate_inner_fold_plan(outer_ids, folds, outer_test_ids=[outer_ids[0]])
    for bad in ([], ["x", "x"]):
        with pytest.raises(ValueError, match="outer_test_ids must be nonempty"):
            th.validate_inner_fold_plan(outer_ids, folds, outer_test_ids=bad)
    with pytest.raises(TypeError):
        th.validate_inner_fold_plan(outer_ids, folds)

