"""Stratified MGMT anchor (src/trackc_calibration.py constants and script 201): published numbers, responder mask, flags.
Synthetic worlds only; no result file is read."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_calibration as cal  # noqa: E402

MIX_WORLD = {"world": 0, "log_r_shift": 0.0, "log_m_scale_non": float(np.log(0.1)), "log_f0_shift": 0.0,
             "B0_center": 0.3, "c": 0.1, "mu": 0.0, "nu": 0.0, "q": 0.35}


def _load():
    spec = importlib.util.spec_from_file_location("stage2_201", REPO / "src" / "201_trackc_stratified_anchor.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_published_stratified_targets_are_in_days_and_the_fraction_interval_is_wilson():
    assert cal.MGMT_UNMETHYLATED_PFS_DAYS[0] == pytest.approx(5.1 * 30.44) and cal.MGMT_UNMETHYLATED_PFS_DAYS[1] == pytest.approx(6.1 * 30.44)
    assert cal.MGMT_METHYLATED_PFS_DAYS[0] == pytest.approx(6.6 * 30.44) and cal.MGMT_METHYLATED_PFS_DAYS[1] == pytest.approx(11.2 * 30.44)
    lo, hi = cal.MGMT_METHYLATED_FRACTION_CI
    assert cal.MGMT_COUNTS == (245, 762) and lo < 245 / 762 < hi
    assert lo == pytest.approx(0.2893, abs=0.001) and hi == pytest.approx(0.3555, abs=0.001)


def test_wilson_interval_matches_known_values():
    lo, hi = cal.wilson_interval(50, 100)
    assert lo == pytest.approx(0.4038, abs=0.001) and hi == pytest.approx(0.5962, abs=0.001)
    lo, hi = cal.wilson_interval(0, 20)
    assert lo == pytest.approx(0.0, abs=1e-9) and 0.1 < hi < 0.2


def test_responder_mask_recovers_the_drawn_classes():
    scen = cal.draw_patients(MIX_WORLD, 4000, [3, 4])
    mask = cal.responder_mask(scen, MIX_WORLD)
    assert abs(mask.mean() - 0.35) < 0.03
    assert scen["m"][mask].min() >= 6.0 - 1e-9 and scen["m"][~mask].max() <= 20 * 0.1 + 1e-9
    with pytest.raises(ValueError):
        cal.responder_mask(scen, {k: v for k, v in MIX_WORLD.items() if k != "q"})


def test_201_evaluate_world_flags_follow_the_published_intervals():
    s = _load()
    cfg = cal.get_config("mixture")
    row = s.evaluate_world(MIX_WORLD, 300, cfg)
    assert 0.25 < row["responder_share_sample"] < 0.45
    lo, hi = cal.MGMT_METHYLATED_FRACTION_CI
    assert row["q_ok"] == (lo <= 0.35 <= hi) and row["q_ok"] is True
    assert row["unmeth_ok"] == (cal.MGMT_UNMETHYLATED_PFS_DAYS[0] <= row["median_non_responders"] <= cal.MGMT_UNMETHYLATED_PFS_DAYS[1])
    assert row["all_ok"] == (row["unmeth_ok"] and row["meth_ok"] and row["q_ok"])
    # responders keep the strong kill range, so they outlast the non-responders in the same trial
    assert row["median_responders"] > row["median_non_responders"] or not np.isfinite(row["median_responders"])


def test_201_analyse_counts_flags_and_lists_strictly_anchored_worlds():
    s = _load()

    def r(world, u, m, q, med_r=250.0, med_n=170.0):
        return {"world": world, "q": q, "median_responders": med_r, "median_non_responders": med_n, "unmeth_ok": u, "meth_ok": m,
                "q_ok": q_ok(q), "all_ok": u and m and q_ok(q)}

    q_ok = lambda q: cal.MGMT_METHYLATED_FRACTION_CI[0] <= q <= cal.MGMT_METHYLATED_FRACTION_CI[1]  # noqa: E731
    rows = [r(0, True, True, 0.32), r(1, True, True, 0.5), r(2, True, False, 0.31), r(3, False, False, 0.33)]
    out = s.analyse(rows, controlled_ids=[0, 2])
    assert out["n_kept_worlds"] == 4 and out["kept"]["n_all_ok"] == 1 and out["kept"]["n_unmeth_ok"] == 3
    assert out["kept"]["n_meth_ok"] == 2 and out["kept"]["n_q_ok"] == 3 and out["kept"]["n_unmeth_and_meth_ok"] == 2
    assert out["controlled"]["n_worlds"] == 2 and out["controlled"]["n_all_ok"] == 1
    assert out["strictly_anchored_world_ids"] == [0] and out["strictly_anchored_and_controlled_world_ids"] == [0]
    assert s.analyse([], [])["n_kept_worlds"] == 0
