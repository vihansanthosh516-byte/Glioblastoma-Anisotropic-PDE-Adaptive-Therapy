"""Cheap dosing rules for the Stage 2 regime map (SIMULATION only). Each rule reads only what a scan shows.

They follow the policy interface of src/trackc_sim.py: reset(scen, rng), act(j, y) -> bool array, advance(j, treat),
and a class attribute uses_filter. NaN in y means no scan this cycle; every rule then keeps its last decision.

    Schedule          open loop: a fixed 0/1 schedule (shared, or one row per patient)
    NadirContainment  induction courses, then treat while the scan shows the burden at least a_on times the smallest
                      earlier scan, and stop once it is back to a_off times that value or lower. The nadir is the
                      smallest noisy reading so far, so it is biased low and the rule fires early when readings are noisy.
"""
from __future__ import annotations

import numpy as np


class Schedule:
    """Treat in the cycles marked in `mask`: shape (N_CYCLES,) shared by all patients, or (N, N_CYCLES) per patient."""
    uses_filter = False

    def __init__(self, mask):
        self.mask = np.asarray(mask, bool)
        if self.mask.ndim not in (1, 2):
            raise ValueError("mask must be 1-D or 2-D")

    def reset(self, scen, rng):
        pass

    def act(self, j, y):
        col = self.mask[j] if self.mask.ndim == 1 else self.mask[:, j]
        return np.broadcast_to(col, y.shape).copy()

    def advance(self, j, treat):
        pass


class NadirContainment:
    uses_filter = False

    def __init__(self, induction, a_on, a_off):
        if not (a_off <= a_on):
            raise ValueError("a_off must not exceed a_on")
        self.induction, self.log_on, self.log_off = int(induction), float(np.log(a_on)), float(np.log(a_off))

    def reset(self, scen, rng):
        self.nadir, self.on = None, None

    def act(self, j, y):
        seen = ~np.isnan(y)
        if j == 0 or self.nadir is None:
            self.nadir = np.where(seen, y, np.inf)
            self.on = np.zeros(y.shape, bool)
        else:
            self.nadir = np.where(seen, np.minimum(self.nadir, y), self.nadir)
        gap = np.where(seen & np.isfinite(self.nadir), y - self.nadir, 0.0)
        new_on = np.where(self.on, gap > self.log_off, gap >= self.log_on)
        self.on = np.where(seen, new_on, self.on)
        out = self.on.copy()
        if j < self.induction:
            out[:] = True
        return out

    def advance(self, j, treat):
        pass
