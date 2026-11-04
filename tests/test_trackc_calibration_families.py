"""Stage 2 run configurations and the mixture world family (src/trackc_calibration.py). Synthetic parameters only."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_calibration as cal  # noqa: E402
import trackc_sim as tc  # noqa: E402

MIX_WORLD = {"world": 0, "log_r_shift": 0.0, "log_m_scale_non": float(np.log(0.1)), "log_f0_shift": 0.0,
             "B0_center": 0.3, "c": 0.1, "mu": 0.0, "nu": 0.0, "q": 0.4}


def test_world_keys_by_family():
    assert cal.world_keys("single") == cal.WORLD_KEYS
    mk = cal.world_keys("mixture")
    assert "q" in mk and "log_m_scale_non" in mk and "log_m_scale" not in mk and mk[-2:] == ["mu", "nu"]
    with pytest.raises(KeyError):
        cal.world_keys("other")


def test_single_family_draws_do_not_change_with_the_family_argument():
    assert cal.draw_worlds(30, 11) == cal.draw_worlds(30, 11, "single")
    assert "q" not in cal.draw_worlds(3, 11)[0]


def test_mixture_worlds_are_stratified_inside_their_axes():
    a = cal.draw_worlds(40, 5, "mixture")
    assert a == cal.draw_worlds(40, 5, "mixture") and a != cal.draw_worlds(40, 6, "mixture")
    for name, lo, hi in cal.MIX_AXES:
        vals = np.array([w[name] for w in a])
        assert (vals >= lo).all() and (vals <= hi).all(), name
        bins = np.floor((vals - lo) / (hi - lo) * 40).astype(int).clip(0, 39)
        assert sorted(bins.tolist()) == list(range(40)), name             # Latin hypercube: one point per stratum
    assert {w["mu"] for w in a} <= set(cal.MU_NU_LEVELS) and {w["nu"] for w in a} <= set(cal.MU_NU_LEVELS)


def test_mixture_world_ranges_use_the_non_responder_scale():
    rr = cal.world_ranges(MIX_WORLD)
    assert rr["m"][0] == pytest.approx(0.6) and rr["m"][1] == pytest.approx(2.0)


def test_mixture_patients_split_into_responders_and_non_responders():
    scen = cal.draw_patients(MIX_WORLD, 4000, [3, 4])
    non_max = tc.RANGES["m"][1] * 0.1                                      # 2.0: top of the non-responder range
    resp = scen["m"] > non_max + 1e-9
    assert abs(resp.mean() - 0.4) < 0.03
    assert scen["m"][resp].min() >= tc.RANGES["m"][0] - 1e-9 and scen["m"][resp].max() <= tc.RANGES["m"][1] + 1e-9
    assert scen["m"][~resp].max() <= non_max + 1e-9


def test_mixture_changes_only_the_kill_strength_of_the_patients():
    w0 = {k: v for k, v in MIX_WORLD.items() if k != "q"}
    base = cal.draw_patients(w0, 500, [9, 9])
    mix = cal.draw_patients(MIX_WORLD, 500, [9, 9])
    for name in ("r", "f0", "B0", "c", "eps", "sigma", "mu", "nu"):
        np.testing.assert_array_equal(mix[name], base[name])
    assert (mix["m"] != base["m"]).any()
    np.testing.assert_array_equal(cal.draw_patients(dict(MIX_WORLD, q=0.0), 500, [9, 9])["m"], base["m"])
    all_resp = cal.draw_patients(dict(MIX_WORLD, q=1.0), 500, [9, 9])["m"]
    assert all_resp.min() >= tc.RANGES["m"][0] - 1e-9


def test_mixture_draw_is_reproducible_and_accepts_an_int_seed():
    a, b = cal.draw_patients(MIX_WORLD, 100, [5, 6]), cal.draw_patients(MIX_WORLD, 100, [5, 6])
    np.testing.assert_array_equal(a["m"], b["m"])
    assert len(cal.draw_patients(MIX_WORLD, 50, 7)["m"]) == 50


def test_more_responders_lengthen_the_standard_arm_median():
    med = lambda q: cal.pfs_stats(*cal.run_arm(cal.draw_patients(dict(MIX_WORLD, q=q), 300, [1, 2]))[:2])["median_interp"]  # noqa: E731
    assert med(1.0) > med(0.0)


def test_config_seed_blocks_are_distinct_and_inside_the_reserved_range():
    seed_keys = ("world_seed", "patient_seed", "d_seed", "e_eval_seed", "e_sel_seed", "i_seed", "j_seed")
    noise_keys = ("noise_c", "noise_d", "noise_e", "noise_i")
    seeds = [c[k] for c in cal.CONFIGS.values() for k in seed_keys]
    noises = [c[k] for c in cal.CONFIGS.values() for k in noise_keys]
    assert len(set(seeds)) == len(seeds) and len(set(noises)) == len(noises)
    assert all(20261801 <= s <= 20262099 for s in seeds) and all(1_810_000 <= n <= 1_899_999 for n in noises)
    for c in cal.CONFIGS.values():                                         # the pairs the manifests register as ranges
        assert c["patient_seed"] == c["world_seed"] + 1 and c["e_sel_seed"] == c["e_eval_seed"] + 1
        assert c["family"] in cal.FAMILY_AXES
    s = cal.CONFIGS["single"]                                              # the run of 2026-11-04 must stay reproducible
    assert (s["world_seed"], s["patient_seed"], s["d_seed"], s["e_eval_seed"], s["e_sel_seed"]) == \
        (20261801, 20261802, 20261803, 20261804, 20261805)
    assert (s["noise_c"], s["noise_d"], s["noise_e"]) == (1_810_000, 1_820_000, 1_830_000)
    assert s["extra_variants"] is False and s["suffix"] == ""


def test_get_config_and_stage_dir():
    cfg = cal.get_config("mixture")
    assert cfg["name"] == "mixture" and cfg["family"] == "mixture"
    assert cal.stage_dir("C", cal.get_config("single")).name == "trackc_stage2_C_calibration"
    assert cal.stage_dir("D", cfg).name == "trackc_stage2_D_negative_controls_mixture"
    assert cal.stage_dir("E", cfg, smoke=True).name == "trackc_stage2_E_regime_map_mixture_smoke"
    with pytest.raises(KeyError):
        cal.get_config("nope")
    cal.get_config("single")["world_seed"] = 1                             # a copy: the table is untouched
    assert cal.CONFIGS["single"]["world_seed"] == 20261801
