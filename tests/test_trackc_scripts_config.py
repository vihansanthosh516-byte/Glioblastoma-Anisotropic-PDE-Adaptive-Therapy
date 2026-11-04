"""Configuration plumbing of the Stage 2 scripts 192 and 193 (run in process on a tiny synthetic mixture family).
The scripts have digit-leading names, so they are loaded by path. Nothing is written outside tmp_path."""
from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_calibration as cal  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location("stage2_" + name, REPO / "src" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _weak_mixture_worlds(n):
    base = {"log_r_shift": 0.3, "log_m_scale_non": float(np.log(0.05)), "log_f0_shift": 0.5, "B0_center": 0.3,
            "c": 0.1, "mu": 0.0, "nu": 0.0, "q": 0.2}
    return [dict(base, world=i, log_r_shift=0.1 * i, q=0.1 + 0.05 * i) for i in range(n)]


def test_192_variants_depend_on_the_configuration():
    s = _load("192_trackc_calibration")
    single = s.variants_for(cal.get_config("single"))
    assert list(single) == list(s.VARIANTS) and "courses3" not in single
    rep = s.variants_for(cal.get_config("single_rep"))
    assert rep["courses3"]["courses"] == 3 and rep["courses4"]["courses"] == 4 and "base" in rep
    assert s.variants_for(cal.get_config("mixture")) == rep


def test_192_job_analyse_and_csv_run_on_the_mixture_family(tmp_path):
    s = _load("192_trackc_calibration")
    cfg = cal.get_config("mixture")
    variants, keys = s.variants_for(cfg), cal.world_keys("mixture")
    worlds = _weak_mixture_worlds(6)
    rows = s._job((0, worlds, 60, cfg, variants))
    assert len(rows) == 6 and set(variants) <= set(rows[0])
    assert rows[0]["courses3"]["mean_courses"] <= 3.0 and rows[0]["courses4"]["mean_courses"] <= 4.0
    res = s.analyse(rows, None, variants, keys)
    assert "courses3" in res["variants"] and "q" in res["axes"] and "log_m_scale" not in res["axes"]
    s.write_csvs(rows, tmp_path, variants, keys)
    header = next(csv.reader(open(tmp_path / "worlds.csv", encoding="utf-8")))
    assert "q" in header and "courses3_pass" in header and "courses4_median_interp" in header
    assert (tmp_path / "kept_worlds.csv").exists()


def test_192_default_arguments_still_describe_the_single_family(tmp_path):
    s = _load("192_trackc_calibration")
    cfg = cal.get_config("single")
    worlds = cal.draw_worlds(6, 3)
    rows = s._job((0, worlds, 40, cfg, s.variants_for(cfg)))
    res = s.analyse(rows)                                                   # no variants, no keys: the old call
    assert set(res["variants"]) == set(s.VARIANTS) and set(res["axes"]) == set(cal.WORLD_KEYS)
    s.write_csvs(rows, tmp_path)
    assert "log_m_scale" in next(csv.reader(open(tmp_path / "worlds.csv", encoding="utf-8")))


def test_193_evaluate_world_and_csv_use_the_family_keys(tmp_path):
    s = _load("193_trackc_negative_controls")
    cfg = cal.get_config("mixture")
    keys = cal.world_keys("mixture")
    w = dict(_weak_mixture_worlds(1)[0], base_median_interp=160.0)
    row = s.evaluate_world(w, 200, cfg, keys)
    assert all(k in row for k in keys) and "dd2_hr" in row and "lm_hr" in row and "controlled" in row
    row["controlled"] = True
    s.write_csvs([row], tmp_path, keys)
    header = next(csv.reader(open(tmp_path / "controlled_worlds.csv", encoding="utf-8")))
    assert header[:1] == ["world"] and "q" in header and "log_m_scale_non" in header
    res = s.analyse([row], keys)
    assert "q" in res["axes"] and res["n_controlled"] == 1


def test_194_analyse_accepts_family_keys():
    s = _load("194_trackc_regime_map")
    worlds = [dict(w, world=i) for i, w in enumerate(_weak_mixture_worlds(5))]
    rows = []
    for w in worlds:
        for ep in s.ENDPOINTS:
            for arm in ("soc6", "front3", "best_fixed", "informed_family_max"):
                for tau in (365, 730):
                    d = 0.0 if arm == "soc6" else 1.0 + w["q"]
                    rows.append({"world": w["world"], "endpoint": ep, "arm": arm, "tau": tau, "rmst": 100.0,
                                 "diff_vs_soc6": d, "ci_lo": d - 1, "ci_hi": d + 1, "regret_cvar10": -d,
                                 "frac_worse": 0.1, "mean_courses": 3.0, "courses_diff_vs_soc6": 0.0 if arm == "soc6" else -3.0,
                                 "frac_beats_informed_730": 0.0, "note": ""})
    res = s.analyse(rows, worlds, cal.world_keys("mixture"))
    block = res["endpoints"]["3b_rano"]["365"]["what_decides_the_gain"]
    assert any("q" in v["spearman"] for v in block.values())
