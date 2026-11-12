"""Forecasts on the amended log measurement scale."""

from math import exp, expm1

TAU = 56.0


def predict_b1(y1, y2):
    return y2


def predict_b2(t1, y1, t2, y2, t3):
    return y2 + (y2 - y1) / (t2 - t1) * (t3 - t2)


def predict_b3(t1, y1, t2, y2, t3):
    return y2 + 0.5 * (y2 - y1) / (t2 - t1) * (t3 - t2)


def fit_m(t1, y1, t2, y2):
    elapsed = t2 - t1
    if y2 < y1:
        return 0.0, (y2 - y1) / (TAU * expm1(-elapsed / TAU))
    return (y2 - y1) / elapsed, 0.0


def predict_m(t1, y1, t2, y2, t3):
    g, k = fit_m(t1, y1, t2, y2)
    elapsed = t3 - t2
    if k == 0:
        return predict_b2(t1, y1, t2, y2, t3)

    return y2 + g * elapsed + k * TAU * exp(-(t2 - t1) / TAU) * expm1(-elapsed / TAU)


def score(pred, obs):
    return abs(pred - obs)

