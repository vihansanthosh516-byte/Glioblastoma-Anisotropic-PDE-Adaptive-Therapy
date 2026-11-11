"""Named particle-filter priors for the integrated Track C benchmark.

See paper/trackc_amendment_pf_prior_support_2026-11-10.md (revision 3).
V3 covers the integrated generator support for all six filter parameters.
NARROW_V2 is the named sensitivity: identical to V3 except r and B0, which keep
the earlier integrated bounds. It equals 208's INTEGRATED_BELIEF_RANGES and is
NOT legacy trackc_sim.RANGES. Legacy defaults elsewhere are unchanged.
"""
import math

F0_LOW = 0.001 * math.exp(-0.7)

PF_PRIOR_INTEGRATED_V3 = {
    "r": (0.00049, 0.124, "log"),
    "m": (0.1, 20.0, "log"),
    "f0": (F0_LOW, 0.9, "log"),
    "c": (0.0, 0.25, "lin"),
    "eps": (0.0, 0.2, "lin"),
    "B0": (0.02, 0.65, "lin"),
}

PF_PRIOR_NARROW_V2 = dict(PF_PRIOR_INTEGRATED_V3, r=(0.001, 0.05, "log"), B0=(0.05, 0.5, "lin"))

__all__ = ["PF_PRIOR_INTEGRATED_V3", "PF_PRIOR_NARROW_V2"]
