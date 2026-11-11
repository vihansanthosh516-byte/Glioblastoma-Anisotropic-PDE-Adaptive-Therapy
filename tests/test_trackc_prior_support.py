"""Outcome-blind check: does the particle-filter prior cover the integrated generator support?

The science-corrected plan (proposed, 2026-11-09, step C2.6) states that the
primary particle-filter prior must cover the full data-generating support.
This test computes the analytic generator support from the sampler's declared
world axes and compares it with the prior in the config and in ``trackc_sim``.
Known gaps are pinned as strict expected failures, like C9-C11 in
``paper/gap_status_2026-11-09.md``. No policy outcome is generated.
Changing the prior is a protocol decision, not a code fix.
"""
import json
import math
from pathlib import Path

import numpy as np
import pytest

import trackc_integrated_sampler as smp
import trackc_sim as tc

ROOT = Path(__file__).resolve().parent.parent
CONFIG = json.loads((ROOT / "configs" / "trackABC_integrated_v2.json").read_text(encoding="utf-8"))
PRIOR = CONFIG["particle_filter_prior"]


def _axis(family, name):
    for axis, lo, hi in smp.WORLD_AXIS_BOUNDS[family]:
        if axis == name:
            return lo, hi
    raise KeyError(name)


def generator_support():
    lo_r, hi_r = _axis("single", "log_r_shift")
    lo_c, hi_c = _axis("single", "B0_center")
    return {
        "r": (tc.RANGES["r"][0] * math.exp(lo_r), tc.RANGES["r"][1] * math.exp(hi_r)),
        "B0": (lo_c - smp.B0_CENTER_HALF_WIDTH, hi_c + smp.B0_CENTER_HALF_WIDTH),
        "m": (tc.RANGES["m"][0] * math.exp(_axis("mixture", "log_m_scale_non")[0]), tc.RANGES["m"][1]),
        "c": tc.RANGES["c"][:2],
        "eps": tc.RANGES["eps"][:2],
    }


def test_config_prior_matches_sim_ranges():
    assert tuple(PRIOR["r_day_minus1_log_uniform"]) == tc.RANGES["r"][:2]
    assert tuple(PRIOR["B0_uniform"]) == tc.RANGES["B0"][:2]


def test_sampled_patients_stay_inside_analytic_support():
    sup = generator_support()
    for family in ("single", "mixture"):
        for seed in range(40):
            world = smp.draw_world_parameters(family, 9_000_000 + seed)
            cohort = smp.draw_patient_cohort(world, 12, 9_100_000 + seed)
            for key in ("r", "B0", "m"):
                assert sup[key][0] - 1e-12 <= cohort[key].min()
                assert cohort[key].max() <= sup[key][1] + 1e-12


@pytest.mark.parametrize("key, prior_key", [("m", "m_log_uniform"), ("c", "c_uniform"), ("eps", "epsilon_uniform")])
def test_prior_covers_support(key, prior_key):
    lo, hi = generator_support()[key]
    assert PRIOR[prior_key][0] <= lo and hi <= PRIOR[prior_key][1]


@pytest.mark.xfail(strict=True, reason="Known gap: PF prior r [0.001, 0.05] vs generator ~[0.0005, 0.123]")
def test_prior_covers_r_support():
    lo, hi = generator_support()["r"]
    assert PRIOR["r_day_minus1_log_uniform"][0] <= lo and hi <= PRIOR["r_day_minus1_log_uniform"][1]


@pytest.mark.xfail(strict=True, reason="Known gap: PF prior B0 [0.05, 0.5] vs generator [0.02, 0.65]")
def test_prior_covers_B0_support():
    lo, hi = generator_support()["B0"]
    assert PRIOR["B0_uniform"][0] <= lo and hi <= PRIOR["B0_uniform"][1]


def test_fraction_outside_prior_is_recorded():
    """Outcome-blind size of the gap on synthetic draws (no policy run)."""
    out = {"r": 0, "B0": 0}
    total = 0
    for family in ("single", "mixture"):
        for seed in range(100):
            world = smp.draw_world_parameters(family, 9_200_000 + seed)
            cohort = smp.draw_patient_cohort(world, 12, 9_300_000 + seed)
            total += 12
            for key, pk in (("r", "r_day_minus1_log_uniform"), ("B0", "B0_uniform")):
                v = cohort[key]
                out[key] += int(np.sum((v < PRIOR[pk][0]) | (v > PRIOR[pk][1])))
    assert out["r"] > 0 and out["B0"] > 0
    print(f"\nOUTSIDE_PRIOR total={total} r={out['r']} B0={out['B0']}")
