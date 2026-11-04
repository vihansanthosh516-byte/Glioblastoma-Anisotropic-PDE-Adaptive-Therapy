"""Script 200 (analysis of the saved outputs of tests C and D): counts, flips, hazard-ratio shares and rank tests on
small synthetic tables. No simulation runs."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_calibration as cal  # noqa: E402


def _load():
    spec = importlib.util.spec_from_file_location("stage2_200", REPO / "src" / "200_trackc_cd_analysis.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _c_row(i, scale, interp, step, other_pass=1):
    lo, hi = cal.TARGET_PFS_DAYS
    return {"world": i, "log_r_shift": 0.0, "log_m_scale": float(np.log(scale)), "log_f0_shift": 0.0, "B0_center": 0.3,
            "c": 0.1, "mu": 0.0, "nu": 0.0, "base_median_interp": interp, "base_median_step": step,
            "base_pass": int(lo <= interp <= hi), "other_pass": other_pass}


def test_analyse_c_compares_the_interpolated_and_the_step_rule():
    s = _load()
    rows = [_c_row(0, 0.1, 150.0, 168.0), _c_row(1, 0.1, 190.0, 168.0), _c_row(2, 0.8, 160.0, 252.0),
            _c_row(3, 0.9, 300.0, 252.0), _c_row(4, 0.05, 100.0, 84.0, other_pass=0)]
    out = s.analyse_c(rows, cal.world_keys("single"))
    assert out["n_worlds"] == 5 and out["interpolated_rule"]["n_pass"] == 2
    st = out["step_rule"]
    assert st["n_pass"] == 2 and st["n_pass_in_both_rules"] == 1 and st["jaccard_with_interpolated"] == pytest.approx(1 / 3)
    assert out["step_median_of_interpolated_passers"] == {"168": 1, "252": 1}
    big = out["worlds_with_scale_ge_line"]
    assert big["n"] == 2 and big["n_pass_interpolated"] == 1 and big["median_interp_of_passers"] == [160.0]
    assert out["variants"]["other"]["n_pass"] == 4 and out["variants"]["base"]["n_pass"] == 2
    assert out["variants"]["base"]["jaccard_with_base"] == 1.0
    assert out["axis_medians"]["kill_scale"]["interpolated_passers"] == pytest.approx(np.median([0.1, 0.8]))
    assert out["kill_scale_axis"] == "log_m_scale"


def test_analyse_c_uses_the_non_responder_scale_in_the_mixture_family():
    s = _load()
    keys = cal.world_keys("mixture")
    rows = []
    for i in range(4):
        r = _c_row(i, 0.1, 150.0 if i < 2 else 300.0, 168.0)
        r["log_m_scale_non"] = r.pop("log_m_scale")
        r["q"] = 0.2 + 0.1 * i
        rows.append(r)
    out = s.analyse_c(rows, keys)
    assert out["kill_scale_axis"] == "log_m_scale_non" and out["interpolated_rule"]["n_pass"] == 2
    assert out["variants"]["base"]["median_q_of_passers"] == pytest.approx(0.25)


def _d_row(i, lm_ok, r_shift, c_med=160.0, re_med=160.0, dd_ok=1, lm_hr=0.9):
    lo, hi = cal.TARGET_PFS_DAYS
    return {"world": i, "log_r_shift": r_shift, "log_m_scale": -2.0, "log_f0_shift": 0.5, "B0_center": 0.3, "c": 0.1,
            "mu": 0.01, "nu": 0.01, "c_median_interp": c_med, "re_median_interp": re_med,
            "re_pass": int(lo <= re_med <= hi), "dd1.5_ok": 1, "dd2_ok": dd_ok, "dd3_ok": 0, "dd2_hr": 0.8,
            "lm_hr": lm_hr, "lm_ok": lm_ok, "controlled": int(lm_ok and dd_ok and lo <= re_med <= hi)}


def test_analyse_d_counts_flips_blumenthal_and_rank_tests():
    s = _load()
    keys = cal.world_keys("single")
    rows = [_d_row(i, 1, -0.5 + 0.01 * i, lm_hr=0.9) for i in range(4)]
    rows += [_d_row(4 + i, 0, 0.5 + 0.01 * i, lm_hr=0.5) for i in range(3)]
    rows += [_d_row(7, 0, 0.6, c_med=144.0, re_med=140.0, lm_hr=0.97)]            # flips: first draw 0.93 d above the edge
    out = s.analyse_d(rows, keys)
    assert out["n_kept_worlds"] == 8 and out["n_still_in_target"] == 7 and out["n_lm_ok"] == 4
    assert len(out["flips"]) == 1 and out["flips"][0]["world"] == 7
    assert out["flips"][0]["first_draw_distance_to_nearest_edge_days"] == pytest.approx(144.0 - cal.TARGET_PFS_DAYS[0])
    assert out["n_controlled"] == 4
    assert out["second_landmark_interval"]["n_worlds_inside"] == 5                   # 0.9 (x4) and 0.97 are inside 0.65-0.98
    rt = out["rank_tests"]["lm_ok"]["log_r_shift"]
    assert rt["auc_flag_higher"] == 0.0 and rt["p_two_sided"] < 0.05 and rt["n_flag"] == 4 and rt["n_not_flag"] == 4
    assert rt["p_bonferroni_over_axes"] >= rt["p_two_sided"]
    assert rt["mean_flag"] < rt["mean_not_flag"]
    assert out["rank_tests"]["lm_ok"]["log_m_scale"] is None                         # no spread: no test
    assert out["controlled_axis_ranges"]["log_r_shift"] == pytest.approx([-0.5, -0.47])
    assert out["controlled_mu_nu_cells"] == {"mu=0.01,nu=0.01": 4}
    assert out["landmark_hr_q05_q25_q50_q75_q95"][2] == pytest.approx(0.9)


def test_analyse_d_handles_an_empty_table():
    assert _load().analyse_d([], cal.world_keys("single")) == {"n_kept_worlds": 0, "primary_dose_column": "dd2"}
