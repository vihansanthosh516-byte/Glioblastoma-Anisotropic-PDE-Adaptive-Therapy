"""Delivered-exposure option of the Track C calibration code (std_caps, array course caps, landmark and dose-dense
controls with a mixed cap, configs single_mix and mixture_mix). Synthetic patients only; no result file is read."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_calibration as cal  # noqa: E402
import trackc_seed_registry as reg  # noqa: E402

WORLD = {"world": 0, "log_r_shift": 0.2, "log_m_scale": float(np.log(0.3)), "log_f0_shift": 0.3, "B0_center": 0.3,
         "c": 0.1, "mu": 0.0, "nu": 0.0}
MIX_WORLD = {"world": 0, "log_r_shift": 0.2, "log_m_scale_non": float(np.log(0.1)), "log_f0_shift": 0.3,
             "B0_center": 0.3, "c": 0.1, "mu": 0.0, "nu": 0.0, "q": 0.35}


def _load(name):
    spec = importlib.util.spec_from_file_location("stage2_" + name, REPO / "src" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_std_caps_int_and_mix36():
    assert cal.std_caps(6, 50) == 6 and cal.std_caps(3, 50) == 3 and np.ndim(cal.std_caps(6, 50)) == 0
    for n in (1, 37, 100, 250, 400, 1000):
        caps = cal.std_caps("mix36", n)
        assert caps.shape == (n,) and set(np.unique(caps)) <= {3, 6}
        assert int((caps == 6).sum()) == (n * 37) // 100            # the exact share for every n, not only for multiples of 100
    caps = cal.std_caps("mix36", 100)
    assert int((caps == 6).sum()) == 37
    sixes = np.flatnonzero(caps == 6)
    assert np.diff(sixes).max() <= 3 and sixes[0] <= 2             # spread evenly over the index, no block of one cap
    # the split depends only on the patient index, so a longer batch starts with the same caps
    assert np.array_equal(cal.std_caps("mix36", 100), cal.std_caps("mix36", 1000)[:100])
    arr = np.array([3, 6, 3])                                       # per-patient caps pass through
    assert np.array_equal(cal.std_caps(arr, 3), arr)
    with pytest.raises(ValueError):
        cal.std_caps(arr, 4)
    with pytest.raises(ValueError):
        cal.std_caps("mix37", 10)


def test_array_caps_equal_the_int_path_when_all_caps_agree():
    scen = cal.draw_patients(WORLD, 80, [7, 1])
    a = cal.run_arm(scen, courses=4, seed=[1, 2])
    b = cal.run_arm(scen, courses=np.full(80, 4), seed=[1, 2])
    for x, y in zip(a, b):
        assert np.array_equal(x, y)
    with pytest.raises(ValueError):
        cal.run_arm(scen, courses=np.full(79, 4), seed=[1, 2])


def test_mixed_caps_give_each_patient_the_result_of_its_own_cap():
    scen = cal.draw_patients(WORLD, 200, [7, 2])
    caps = cal.std_caps("mix36", 200)
    mix = cal.run_arm(scen, courses=caps, seed=[1, 2])
    three = cal.run_arm(scen, courses=3, seed=[1, 2])
    six = cal.run_arm(scen, courses=6, seed=[1, 2])
    low, high = caps == 3, caps == 6
    assert np.array_equal(mix[0][low], three[0][low]) and np.array_equal(mix[0][high], six[0][high])
    assert np.array_equal(mix[2][low], three[2][low]) and np.array_equal(mix[2][high], six[2][high])
    assert mix[2][low].max() <= 3 and mix[2][high].max() <= 6


def test_landmark_control_default_is_unchanged_and_mix_uses_only_six_course_patients():
    scen = cal.draw_patients(WORLD, 400, [7, 3])
    old = cal.landmark_control(scen, seed=[1, 2])                       # six courses for everyone, as before
    same = cal.landmark_control(scen, seed=[1, 2], courses=6)
    assert old["hr"] == same["hr"] and old["n_landmark"] == same["n_landmark"] and old["n_total"] == 400
    mix = cal.landmark_control(scen, seed=[1, 2], courses="mix36")
    caps = cal.std_caps("mix36", 400)
    assert mix["n_total"] == 400 and mix["n_landmark"] <= int((caps == 6).sum())
    assert mix["n_landmark"] < old["n_landmark"]
    none = cal.landmark_control(scen, seed=[1, 2], courses=3)           # nobody completes 6 cycles
    assert none["n_landmark"] == 0 and np.isnan(none["hr"])


def test_dose_dense_control_runs_with_a_mixed_cap():
    scen = cal.draw_patients(WORLD, 300, [7, 4])
    caps = cal.std_caps("mix36", 300)
    out = cal.dose_dense_control(scen, mult=2.0, seed=[1, 2], courses=caps)
    assert np.isfinite(out["hr"]) and out["hr"] < 1.0                   # twice the drug effect cannot lengthen PFS here
    again = cal.dose_dense_control(scen, mult=2.0, seed=[1, 2], ref=cal.run_arm(scen, courses=caps, seed=[1, 2]),
                                   courses=caps)
    assert again["hr"] == out["hr"]


def test_new_configs_are_defined_with_distinct_dev_seeds():
    base = cal.get_config("single")
    assert base["std_courses"] == 6 and cal.get_config("single_rep")["std_courses"] == 6
    assert cal.get_config("mixture")["std_courses"] == 6
    sm, mm = cal.get_config("single_mix"), cal.get_config("mixture_mix")
    assert sm["family"] == "single" and mm["family"] == "mixture" and sm["std_courses"] == mm["std_courses"] == "mix36"
    assert sm["extra_variants"] and mm["extra_variants"] and mm["k_seed"] and "k_seed" not in sm
    seed_keys = ("world_seed", "patient_seed", "d_seed", "e_eval_seed", "e_sel_seed", "i_seed", "j_seed", "k_seed")
    noise_keys = ("noise_c", "noise_d", "noise_e", "noise_i", "noise_k")
    taken_s, taken_n = {}, {}
    for name in cal.CONFIGS:
        cfg = cal.get_config(name)
        for k in seed_keys:
            if k in cfg:
                assert cfg[k] not in taken_s, (name, k, taken_s[cfg[k]])
                taken_s[cfg[k]] = (name, k)
                assert reg.STAGE2_DEV_SCENARIO_BLOCK[0] <= cfg[k] <= reg.STAGE2_DEV_SCENARIO_BLOCK[1]
        for k in noise_keys:
            if k in cfg:
                assert cfg[k] not in taken_n, (name, k, taken_n[cfg[k]])
                taken_n[cfg[k]] = (name, k)
                assert reg.STAGE2_DEV_NOISE_BLOCK[0] <= cfg[k] <= reg.STAGE2_DEV_NOISE_BLOCK[1]
    suffixes = [cal.get_config(n)["suffix"] for n in cal.CONFIGS]
    assert len(set(suffixes)) == len(suffixes)


def test_192_variants_and_job_for_a_mixed_exposure():
    s = _load("192_trackc_calibration")
    cfg = cal.get_config("single_mix")
    variants = s.variants_for(cfg)
    assert {"courses3", "courses4", "courses6"} <= set(variants) and variants["courses6"]["courses"] == 6
    assert "courses6" not in s.variants_for(cal.get_config("single_rep"))
    base = s.base_courses(cfg, 3 * 100, 100)
    assert base.shape == (300,) and int((base == 6).sum()) == 111
    assert s.base_courses(cal.get_config("single"), 300, 100) == 6
    rows = s._job((0, cal.draw_worlds(3, 5), 100, cfg, variants))
    for r in rows:
        # the cap is a planned maximum: the share that actually receives 6 courses cannot exceed the 37% with a cap of 6
        assert 0.0 <= r["base"]["frac_courses_ge6"] <= 0.37 + 1e-9 and 0.0 <= r["base"]["median_courses"] <= 6.0
        assert r["courses6"]["frac_courses_ge6"] >= r["base"]["frac_courses_ge6"]
        assert r["courses3"]["mean_courses"] <= 3.0 and r["courses6"]["mean_courses"] <= 6.0
        assert r["courses3"]["mean_courses"] <= r["base"]["mean_courses"] <= r["courses6"]["mean_courses"] + 1e-9
        assert r["courses3"]["median_interp"] <= r["courses6"]["median_interp"] or not np.isfinite(r["courses6"]["median_interp"])


def test_192_analyse_reports_the_courses_actually_delivered():
    s = _load("192_trackc_calibration")
    cfg = cal.get_config("single_mix")
    variants = s.variants_for(cfg)
    rows = s._job((0, cal.draw_worlds(4, 8), 100, cfg, variants))
    res = s.analyse(rows, None, variants, cal.world_keys("single"))
    d = res["variants"]["base"]["delivered_courses"]
    assert len(d["median_courses_q25_50_75_all_worlds"]) == 3 and max(d["frac_ge6_q25_50_75_all_worlds"]) <= 0.37 + 1e-9
    assert res["variants"]["courses3"]["delivered_courses"]["frac_ge6_q25_50_75_all_worlds"] == [0.0, 0.0, 0.0]


def test_193_evaluate_world_with_a_mixed_exposure():
    s = _load("193_trackc_negative_controls")
    cfg = cal.get_config("mixture_mix")
    w = dict(MIX_WORLD, base_median_interp=160.0)
    row = s.evaluate_world(w, 200, cfg, cal.world_keys("mixture"))
    caps = cal.std_caps("mix36", 200)
    assert row["lm_n"] <= int((caps == 6).sum()) and "dd2_hr" in row and "controlled" in row


def test_194_reference_follows_the_exposure_and_the_new_arm_is_a_valid_reference():
    s = _load("194_trackc_regime_map")
    assert s.reference_for(cal.get_config("single")) == "soc6" == s.reference_for(cal.get_config("mixture"))
    assert s.reference_for(cal.get_config("single_mix")) == "std_capmix" == s.reference_for(cal.get_config("mixture_mix"))
    cfg = cal.get_config("single_mix")
    s.N_EVAL, s.N_SEL = 60, 40
    rows = s.evaluate_world(dict(WORLD, base_median_interp=160.0), 60, 40, cfg)
    ref_rows = [r for r in rows if r["arm"] == "std_capmix"]
    assert ref_rows and all(r["diff_vs_ref"] == 0.0 and r["courses_diff_vs_ref"] == 0.0 for r in ref_rows)
    soc6 = [r for r in rows if r["arm"] == "soc6" and r["endpoint"] == "3b_rano" and r["tau"] == 730][0]
    assert soc6["courses_diff_vs_ref"] > 0                                 # six courses for all is more drug than 3-or-6
    plain = s.evaluate_world(dict(WORLD, base_median_interp=160.0), 60, 40, cal.get_config("single"))
    assert not [r for r in plain if r["arm"] == "std_capmix"]
    assert all(r["diff_vs_ref"] == 0.0 for r in plain if r["arm"] == "soc6")


def test_exposure_protocol_leaves_default_protocols_unchanged():
    assert cal.exposure_protocol(cal.get_config("single")) == {}
    assert cal.exposure_protocol(cal.get_config("single"), always=True) == {"std_courses": 6}
    assert cal.exposure_protocol(cal.get_config("mixture_mix")) == {"std_courses": "mix36",
                                                                    "mix36_six_fraction": cal.MIX36_SIX_FRACTION}


def test_scripts_pass_the_course_cap_explicitly():
    """A call that leaves out `courses` would silently use the default of 6 and ignore the configuration."""
    import ast
    guarded = {"run_arm", "dose_dense_control", "landmark_control", "best_response"}
    missing, checked = [], 0
    for path in sorted((REPO / "src").glob("[12][0-9][0-9]_trackc_*.py")):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.Call):
                f = node.func
                name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", None)
                if name in guarded:
                    checked += 1
                    if not any(k.arg == "courses" for k in node.keywords):
                        missing.append(f"{path.name}:{node.lineno} {name}")
    assert checked >= 10 and missing == [], missing


def test_199_best_response_accepts_per_patient_caps():
    s = _load("199_trackc_response_depth")
    scen = cal.draw_patients(WORLD, 120, [9, 1])
    assert np.array_equal(s.best_response(scen, courses=5), s.best_response(scen, courses=np.full(120, 5)))
    mixed = s.best_response(scen, courses=cal.std_caps("mix36", 120))
    low, high = s.best_response(scen, courses=3), s.best_response(scen, courses=6)
    assert np.all(mixed >= high - 1e-12) and np.all(mixed <= low + 1e-12)
    row = s.evaluate_world(WORLD, 120, 9, "mix36")
    assert 0.0 <= row["frac_best_le_0.5"] <= 1.0


def test_201_evaluate_world_with_a_mixed_exposure():
    s = _load("201_trackc_stratified_anchor")
    cfg = cal.get_config("mixture_mix")
    row = s.evaluate_world(MIX_WORLD, 300, cfg)
    assert 0.2 < row["responder_share_sample"] < 0.5 and "all_ok" in row
    assert row["median_responders"] > row["median_non_responders"] or not np.isfinite(row["median_responders"])
