import numpy as np

from src import fk_two_history_sensitivity_screen as sc
from src import fk_two_history_shape as th

N = 24


def test_truth_regimes_behave():
    S = sc.support(N)
    v = {r: [int(sc.truth_mask(r, t, 3, N, S).sum()) for t in (10, 18, 34)] for r in sc.REGIMES}
    assert v["growth"][0] < v["growth"][2]
    assert v["shrinkage"][0] > v["shrinkage"][2]
    assert abs(v["plateau"][2] - v["plateau"][1]) <= 0.15 * max(v["plateau"][1], 1)


def test_case_structure_separate_fields_and_constants_restored():
    before = (th.SHRINK_TAU_DAYS, th.MAX_FACTOR, th.MAX_SD_DELTA_MM)
    c = sc.run_case("growth", 8, 8, 0, n=N, grid=False, n_anchors=32)
    assert (th.SHRINK_TAU_DAYS, th.MAX_FACTOR, th.MAX_SD_DELTA_MM) == before
    assert set(c["arms"]) == {"copy_last", "last_mask_sd", "control_default"}
    assert c["arms"]["copy_last"]["dice_true_vol_boundary"] is None
    for arm in ("last_mask_sd", "control_default"):
        for f in sc.DICE_FIELDS:
            assert 0.0 <= c["arms"][arm][f] <= 1.0
    assert c["arms"]["last_mask_sd"]["factor"] == 0.0 and c["arms"]["control_default"]["factor"] > 0


def test_deterministic_and_summary():
    a = sc.run_screen(regimes=("growth", "shrinkage"), hist=(8,), fcs=(4,), seeds=(0, 1), n=N, grid=False)
    b = sc.run_screen(regimes=("growth", "shrinkage"), hist=(8,), fcs=(4,), seeds=(0, 1), n=N, grid=False)
    assert a["cases"] == b["cases"]
    e = a["by_regime"]["growth"]
    assert e["n_cases"] == 2
    ref = e["arms"]["last_mask_sd"]["dice_true_vol_pre"]
    assert ref["paired_diff_vs_last_mask_sd_mean"] == 0.0
    assert e["arms"]["copy_last"]["dice_true_vol_boundary"] is None


def test_grid_size():
    assert len(sc.control_settings(True)) == 1 + 3 * 3 * 2
