"""Small survival statistics for simulated Track C outcomes (SIMULATION statistics only).

Model outcomes are times to progression with administrative censoring at the horizon, and many ties
when progression is detected only at scan visits. This module gives:
    km_median     Kaplan-Meier median (np.inf if the curve never reaches 0.5)
    km_median_interp  the same with linear interpolation between event times (for coarse scan grids)
    km_survival   Kaplan-Meier survival probability at a time (events at that time count as events)
    cox_hr        two-group Cox hazard ratio, Efron ties, Wald 95% interval
It exists so a simulated contrast can be put on the same scale as the published hazard ratios.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar


def _prep(time, event):
    t = np.asarray(time, float)
    e = np.asarray(event, bool)
    if t.shape != e.shape:
        raise ValueError("time and event must have the same shape")
    return t, e


def km_survival(time, event, at):
    """Kaplan-Meier survival probability at time `at`."""
    t, e = _prep(time, event)
    s = 1.0
    for u in np.unique(t[e & (t <= at)]):
        d = np.sum(e & (t == u))
        n = np.sum(t >= u)
        s *= 1.0 - d / n
    return float(s)


def km_median(time, event):
    """Smallest time at which the Kaplan-Meier curve is at or below 0.5; np.inf if it never gets there."""
    t, e = _prep(time, event)
    s = 1.0
    for u in np.unique(t[e]):
        d = np.sum(e & (t == u))
        n = np.sum(t >= u)
        s *= 1.0 - d / n
        if s <= 0.5 + 1e-12:
            return float(u)
    return float("inf")


def km_median_interp(time, event):
    """Kaplan-Meier median with linear interpolation between event times (the curve starts at S = 1 at t = 0).
    It treats events as spread evenly inside each gap between scans, so a coarse scan grid does not force the
    median onto a scan day. np.inf if the curve never reaches 0.5."""
    t, e = _prep(time, event)
    s_prev, t_prev = 1.0, 0.0
    for u in np.unique(t[e]):
        d = np.sum(e & (t == u))
        n = np.sum(t >= u)
        s = s_prev * (1.0 - d / n)
        if s <= 0.5 + 1e-12:
            return float(t_prev + (u - t_prev) * (s_prev - 0.5) / (s_prev - s))
        s_prev, t_prev = s, float(u)
    return float("inf")


def _neg_loglik(beta, t, e, g):
    """Negative Efron partial log-likelihood for one binary covariate g (0/1)."""
    theta = np.exp(beta)
    ll = 0.0
    for u in np.unique(t[e]):
        at_risk = t >= u
        n0 = np.sum(at_risk & (g == 0))
        n1 = np.sum(at_risk & (g == 1))
        ev = e & (t == u)
        d0, d1 = np.sum(ev & (g == 0)), np.sum(ev & (g == 1))
        d = d0 + d1
        a_all = n0 + theta * n1
        b_ev = d0 + theta * d1
        ll += beta * d1
        for l in range(int(d)):
            ll -= np.log(a_all - (l / d) * b_ev)
    return -ll


def cox_hr(time, event, group):
    """Hazard ratio of group 1 relative to group 0 (Cox model, Efron ties) with a Wald 95% interval.
    Returns (hr, (low, high)). Needs at least one event in each group."""
    t, e = _prep(time, event)
    g = np.asarray(group, int)
    if not (e & (g == 0)).any() or not (e & (g == 1)).any():
        raise ValueError("each group needs at least one event")
    res = minimize_scalar(_neg_loglik, bounds=(-6.0, 6.0), args=(t, e, g), method="bounded",
                          options={"xatol": 1e-9})
    b = float(res.x)
    h = 1e-4
    f0, fp, fm = _neg_loglik(b, t, e, g), _neg_loglik(b + h, t, e, g), _neg_loglik(b - h, t, e, g)
    info = (fp - 2.0 * f0 + fm) / h ** 2
    se = float(np.sqrt(1.0 / info)) if info > 0 else float("nan")
    return float(np.exp(b)), (float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se)))
