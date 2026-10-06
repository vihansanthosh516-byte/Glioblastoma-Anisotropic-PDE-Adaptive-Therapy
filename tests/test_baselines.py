"""Unit tests for the baseline helpers in src/98_baseline_ladder.py (no patient data needed)."""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
spec = importlib.util.spec_from_file_location("baseline_ladder", REPO / "src" / "98_baseline_ladder.py")
bl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bl)


def ball(n=40, r=6):
    z, y, x = np.meshgrid(*[np.arange(n)] * 3, indexing="ij")
    return ((x - n // 2) ** 2 + (y - n // 2) ** 2 + (z - n // 2) ** 2) < r * r


def test_reshape_hits_target_volume_and_nests():
    S = ball()
    brain = np.ones_like(S)
    sd = bl.signed_domain_distance(S, brain)
    n0 = int(S.sum())
    big = bl.reshape_to_volume(S, sd, n0 * 2)
    small = bl.reshape_to_volume(S, sd, n0 // 2)
    assert abs(int(big.sum()) - 2 * n0) <= 0.05 * n0 * 2
    assert abs(int(small.sum()) - n0 // 2) <= 0.05 * n0
    assert (small & ~S).sum() == 0 and (S & ~big).sum() == 0          # shrink stays inside, growth contains S
    same = bl.reshape_to_volume(S, sd, n0)
    assert (same == S).mean() > 0.99


def test_reshape_stays_in_brain_domain():
    S = ball()
    brain = np.zeros_like(S)
    brain[:, :, :22] = True
    brain |= S
    sd = bl.signed_domain_distance(S, brain)
    g = bl.reshape_to_volume(S, sd, int(S.sum()) * 3)
    assert (g & ~brain).sum() == 0


def test_empty_cases():
    S = np.zeros((10, 10, 10), bool)
    sd = bl.signed_domain_distance(S, np.ones_like(S))
    assert not bl.reshape_to_volume(S, sd, 50).any()
    S = ball(20, 4)
    sd = bl.signed_domain_distance(S, np.ones_like(S))
    assert not bl.reshape_to_volume(S, sd, 0).any()


def test_gompertz_recovers_known_curve_and_needs_three_scans():
    K, V0, a = 5000.0, 100.0, 0.02
    days = [0.0, 30.0, 60.0]
    vols = [np.exp(math.log(K + 1) + (math.log(V0 + 1) - math.log(K + 1)) * math.exp(-a * d)) - 1 for d in days]
    pred = bl.gompertz_forecast(days, vols, K, 90.0)
    true = math.exp(math.log(K + 1) + (math.log(V0 + 1) - math.log(K + 1)) * math.exp(-a * 90.0)) - 1
    assert abs(pred / true - 1) < 0.01
    assert bl.gompertz_forecast(days[:2], vols[:2], K, 90.0) is None
    assert bl.gompertz_forecast([0, 30, 60], [0, 10, 20], K, 90.0) is None   # empty first scan
