"""Synthetic tests for the pilot design rule. Made-up numbers; no policy outcome."""
import numpy as np

import trackc_dev_pilot_design as dp

SCREEN = {"rows": [{"family": f, "sd_days": s, "n_required": n}
                   for f in ("normal", "skewed") for s, n in ((20, 50), (40, 200), (80, 800))]}


def test_variance_components_recover_known_values():
    rng = np.random.default_rng(0)
    tau, sigma, m = 10.0, 40.0, 12
    worlds = [rng.normal(rng.normal(0, tau), sigma, m) for _ in range(4000)]
    t2, s2, m0 = dp.variance_components(worlds)
    assert m0 == m
    assert abs(s2 - sigma**2) / sigma**2 < 0.03
    assert abs(t2 - tau**2) / tau**2 < 0.25


def test_n_for_sd_interpolates_and_refuses_above_table():
    assert dp.n_for_sd(30.0, SCREEN) == 125
    assert dp.n_for_sd(90.0, SCREEN) is None


def _recs(tau, sigma, m=12, k=20, secs=100.0):
    rng = np.random.default_rng(1)
    return [{"family": fam, "phase": "dev_pilot", "patients_per_cohort": m, "elapsed_seconds": secs,
             "contrasts": {"minimax": {"per_patient": rng.normal(rng.normal(0, tau), sigma, m).tolist()}}}
            for fam in ("single", "mixture") for _ in range(k)]


def test_patient_noise_dominant_prefers_larger_cohorts_than_world_noise_dominant():
    timing = (300.0, 10.0)  # large fixed cost per world, small cost per patient
    noisy_patients = dp.design(_recs(tau=2.0, sigma=60.0), SCREEN, timing=timing)
    noisy_worlds = dp.design(_recs(tau=30.0, sigma=5.0), SCREEN, timing=timing)
    assert noisy_patients["choice_m"] is not None and noisy_worlds["choice_m"] is not None
    assert noisy_patients["choice_m"] > noisy_worlds["choice_m"]


def test_proportional_time_picks_smallest_cohort_and_says_so():
    out = dp.design(_recs(tau=2.0, sigma=60.0), SCREEN)
    assert out["choice_m"] == min(dp.M_GRID)
    assert out["time_model_source"].startswith("assumed")
