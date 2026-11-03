"""Synthetic tests for the protocol-complete LUMIERE wrapper. No patient data."""
from dataclasses import replace

import numpy as np

from src import fk_lumiere_protocol as lp
from tests.test_fk_lumiere_runner import FOLDS, SHAPE, TINY, FakeProvider, _cases, _folds


def _spread_cases():
    # forecast gap independent of fold (fold = i % 5, gap = 11 + (i // 5) % 5)
    return [replace(c, forecast_weeks=11 + (i // 5) % 5) for i, c in enumerate(_cases())]


def _factory(log):
    def make(rho_b, d_b):
        log.append((tuple(rho_b), tuple(d_b)))
        return FakeProvider()
    return make


def test_three_prior_runs_share_patients_folds_and_film_cache_and_gap_refit_is_separate():
    log = []
    cases = _spread_cases()
    out = lp.run_lumiere_protocol(cases, _folds(cases), provider_factory=_factory(log), recipe=TINY, n_bootstrap=200)
    assert out["status"] == "completed (synthetic)"
    assert log[:3] == [tuple(map(tuple, v)) for v in lp.PRIOR_SCENARIOS.values()] and log[3] == log[0]
    runs = out["runs"]
    ids = set(runs["primary"]["dice"]["candidate"])
    assert all(set(r["dice"]["candidate"]) == ids for r in runs.values())
    for name, r in runs.items():
        for f, info in r["folds"].items():
            assert info["test_ids"] == runs["primary"]["folds"][f]["test_ids"]
            assert not set(info["training_ids"]) & set(info["test_ids"])
        reused = [i["film_reused"] for i in r["folds"].values()]
        assert all(reused) if name != "primary" else not any(reused)
    gap = out["central_gap_run"]
    central = {c.patient_id for c in cases if c.forecast_weeks in lp.CENTRAL_GAP_WEEKS}
    assert set(gap["dice"]["candidate"]) == central
    assert not any(i["film_reused"] for i in gap["folds"].values())  # own cache
    for info in gap["folds"].values():
        assert set(info["training_ids"]) <= central  # edge-gap patients never fit
    assert set(out["central_gap_excluded_ids"]) == {c.patient_id for c in cases} - central
    for f, info in gap["folds"].items():
        assert set(info["test_ids"]) == {p for p in central if FOLDS[p] == f}


def test_primary_stop_stops_every_run():
    log = []
    cases = _cases()
    cases = [replace(c, qc_passed=False) if i < 3 else c for i, c in enumerate(cases)]
    out = lp.run_lumiere_protocol(cases, _folds(cases), provider_factory=_factory(log), recipe=TINY)
    assert out["status"].startswith("stopped in primary") and list(out["runs"]) == ["primary"] and len(log) == 1


def test_frozen_settings_match_draft():
    assert lp.CALIBRATION_MINIMUM_PATIENTS == 8 and lp.CALIBRATION_FOLDS == 5
    assert lp.INITIAL_BOUNDARY_COVARIANCE == {"sigma_b_mm": 0.8, "tau_b_mm": 0.3, "ell_mm": 6.0}
    assert (lp.FORECAST_SEED, lp.BOUNDARY_CALIBRATION_SEED) == (20261028, 20261029)
    assert lp.CENTRAL_GAP_WEEKS == (12, 13, 14)
    prov = lp.default_provider_factory(*lp.PRIOR_SCENARIOS["primary"])
    assert prov.seed == 20261028 and prov.calibration_minimum == 8 and prov.settings["n_draws"] == 200
