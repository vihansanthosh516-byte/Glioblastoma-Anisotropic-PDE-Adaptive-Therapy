"""
Independent reference solver for the Track C simulator (audit C1, 2026-10-26).
=============================================================================
Same equations as trackc_sim._deriv, integrated with SciPy solve_ivp (DOP853, tight
tolerances) piecewise over each cycle's on (0-5 d) and off (5-28 d) segments, with a terminal
event at S + R = 1.2 B0. SciPy is BSD-3-Clause. This checks numerics only; it says nothing
about whether the model is biologically right.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

import trackc_sim as tc


def _rhs(p, u):
    def f(_t, y):
        dS, dR = tc._deriv(y[0], y[1], p, u)
        return [dS, dR]
    return f


def reference_ttp(p, schedule, n_cycles=tc.N_CYCLES, rtol=1e-10, atol=1e-13, method="DOP853"):
    """Capped time to S + R >= 1.2 B0 for one scenario (scalar dict p) and a 0/1 schedule.
    Treatment is not given after progression (the event ends the run)."""
    p = {k: float(np.ravel(v)[0]) for k, v in p.items()}
    level = tc.PROG_FACTOR * p["B0"]
    y = np.array([(1 - p["f0"]) * p["B0"], p["f0"] * p["B0"]])

    def event(_t, yy):
        return yy[0] + yy[1] - level
    event.terminal, event.direction = True, 1

    for j in range(n_cycles):
        t0 = j * tc.CYCLE_DAYS
        for a, b, on in ((0.0, tc.COURSE_DAYS, True), (tc.COURSE_DAYS, tc.CYCLE_DAYS, False)):
            u = 1.0 if (on and schedule[j]) else 0.0
            sol = solve_ivp(_rhs(p, u), (t0 + a, t0 + b), y, method=method, rtol=rtol, atol=atol,
                            events=event)
            if sol.t_events[0].size:
                return float(sol.t_events[0][0])
            y = sol.y[:, -1]
    return float(n_cycles * tc.CYCLE_DAYS)
