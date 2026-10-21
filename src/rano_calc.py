# NOTE (2026-10-19/20 audit): current RANO display/calculator; model-based proxy sensitivity span from one centre and two tools, not formal noise bounds; not validated and not for clinical use. The 500 mm3 floor is a project operational choice, and backward mode assumes a flat prior with normal error.
"""Probabilistic RANO calculator (Amendment 17, analysis_plan_v1.md).

Model-based proxy sensitivity span from one centre (LUMIERE, 2 tools); NOT for clinical use. Not validated on patients.
Volumes in mm3. s = SD of the log-ratio (change) error, ln units, as in A15/A16.
Thresholds: progression ratio >= 1.40, partial response ratio <= 0.35 (RANO 2.0). The 500 mm3 floor is an operational
choice of this project, not a RANO 2.0 rule.
Backward mode uses a flat prior on the true log ratio and a normal error: a flat-prior, normal-error approximation.
"""
from __future__ import annotations

import json
import math
from pathlib import Path as _Path

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
A15_JSON = PROJECT_ROOT / "output" / "rano_noise" / "rano_tool_dependence.json"

PD_RATIO, PR_RATIO = 1.40, 0.35
LN_PD, LN_PR = math.log(PD_RATIO), math.log(PR_RATIO)
FLOOR_MM3 = 500.0
LABEL = "model-based proxy sensitivity scenarios from one centre (LUMIERE, 2 tools); not formal noise bounds; not for clinical use"


def _phi(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def classify(v0: float, v1: float) -> str:
    """Deterministic RANO 2.0 volumetric class of the ratio v1/v0 (no noise)."""
    r = v1 / v0
    return "PD" if r >= PD_RATIO else ("PR" if r <= PR_RATIO else "SD")


def _a15():
    return json.loads(A15_JSON.read_text())["TEST2_mdc_vs_progression_line"]


def size_band(v0: float, v1: float) -> int:
    """A15 size band (0, 1, 2) from the geometric mean of the two volumes."""
    e = _a15()["band_edges_mm3"]
    g = math.sqrt(v0 * v1)
    return int(g >= e[0]) + int(g >= e[1])


def bracket(band: int) -> tuple[float, float]:
    """Legacy API: (min, max) of two A15 change-SD proxy estimates, not formal bounds."""
    b = _a15()["bands"][str(band)]
    lo, hi = sorted((b["s_change_upper_bound"], b["s_change_lower_bound"]))
    return float(lo), float(hi)


def backward(v0: float, v1: float, s: float) -> dict:
    """P(true change is PD / PR / between) given the observed volumes (flat prior, normal error)."""
    m = math.log(v1 / v0)
    if s <= 0:
        c = classify(v0, v1)
        return {"P_PD": float(c == "PD"), "P_PR": float(c == "PR"), "P_SD": float(c == "SD")}
    p_pd = 1.0 - _phi((LN_PD - m) / s)
    p_pr = _phi((LN_PR - m) / s)
    return {"P_PD": p_pd, "P_PR": p_pr, "P_SD": 1.0 - p_pd - p_pr}


def forward(true_ratio: float, s: float) -> dict:
    """P(observed class is PD / PR / SD) given a true volume ratio (observed log ratio ~ N(ln true, s^2))."""
    m = math.log(true_ratio)
    if s <= 0:
        c = "PD" if true_ratio >= PD_RATIO else ("PR" if true_ratio <= PR_RATIO else "SD")
        return {"P_PD": float(c == "PD"), "P_PR": float(c == "PR"), "P_SD": float(c == "SD")}
    p_pd = 1.0 - _phi((LN_PD - m) / s)
    p_pr = _phi((LN_PR - m) / s)
    return {"P_PD": p_pd, "P_PR": p_pr, "P_SD": 1.0 - p_pd - p_pr}


def mdc(s: float) -> dict:
    """Minimal detectable change at 95%: 1.96 s in ln units, and as a percent change."""
    x = 1.96 * s
    return {"mdc95_ln": x, "mdc95_pct": math.expm1(x) * 100.0}


def rano_calc(v0: float, v1: float, s: float | None = None) -> dict:
    """Full output for two volumes. If s is None the A15 bracket for the size band is used (low and high)."""
    below = bool(v0 < FLOOR_MM3 or v1 < FLOOR_MM3)
    out = {"label": LABEL, "v0_mm3": v0, "v1_mm3": v1, "ratio": v1 / v0,
           "below_measurability_floor": below,
           "deterministic_class": "not measurable (below 500 mm3 operational floor)" if below else classify(v0, v1)}
    if s is not None:
        out["s"] = s
        out["backward"] = backward(v0, v1, s)
        out.update(mdc(s))
    else:
        band = size_band(v0, v1)
        lo, hi = bracket(band)
        out["size_band"] = band
        out["s_bracket"] = [lo, hi]
        out["backward_at_low_s"] = backward(v0, v1, lo)
        out["backward_at_high_s"] = backward(v0, v1, hi)
        out["mdc_low"] = mdc(lo)
        out["mdc_high"] = mdc(hi)
    return out


if __name__ == "__main__":
    import sys
    a = [float(x) for x in sys.argv[1:]]
    if len(a) not in (2, 3):
        raise SystemExit("usage: python src/rano_calc.py V0_mm3 V1_mm3 [s]")
    print(json.dumps(rano_calc(a[0], a[1], a[2] if len(a) == 3 else None), indent=2))
