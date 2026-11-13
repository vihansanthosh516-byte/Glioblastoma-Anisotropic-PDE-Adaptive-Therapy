"""Outcome-blind precision screen for the Track C v2 world count (C-6).

For each assumed SD of world-level paired differences, find the smallest number of
worlds per family at which the GE-betting interval (delta 0.025 per family contrast,
support [-790, 790]) has half-width <= 14 days in at least 80% of synthetic repeats.
Synthetic world differences only; no policy outcome is read or produced.

The development pilot (src/trackc_dev_pilot_design.py) measures the real SD and reads
the table this script writes to choose cohort size and world count (decision D3).

  python -I src/trackc_precision_screen.py [--reps 200] [--out output/trackc_precision_screen]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
from trackc_ge_betting import ge_betting_interval  # noqa: E402

SUPPORT = (-790.0, 790.0)
DELTA = 0.025
TARGET = 14.0
QUANTILE = 0.80
SEED = 20261113
SD_GRID = (15.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 120.0)
FAMILIES = ("normal", "skewed")
N_CAP = 20000


def draw(family: str, sd: float, n: int, rng) -> np.ndarray:
    if family == "normal":
        x = rng.normal(10.0, sd, n)
    elif family == "skewed":
        k = 2.0
        x = rng.gamma(k, sd / np.sqrt(k), n) - sd * np.sqrt(k) + 10.0
    else:
        raise ValueError(family)
    return np.clip(x, *SUPPORT)


def half_width_quantile(family, sd, n, reps, rng) -> float:
    widths = []
    for _ in range(reps):
        iv = ge_betting_interval(draw(family, sd, n, rng), support=SUPPORT, delta=DELTA)
        widths.append(np.inf if iv is None else (iv[1] - iv[0]) / 2.0)
    return float(np.quantile(widths, QUANTILE))


def required_n(family, sd, reps, seed) -> dict:
    """Smallest n (to 2% resolution) meeting the rule, by doubling then bisection."""
    def ok(n):
        rng = np.random.Generator(np.random.PCG64DXSM([seed, int(sd * 10), n]))
        return half_width_quantile(family, sd, n, reps, rng) <= TARGET
    lo, hi = 1, 16
    while not ok(hi):
        lo, hi = hi, hi * 2
        if hi > N_CAP:
            return {"family": family, "sd_days": sd, "n_required": None, "note": f"above cap {N_CAP}"}
    while hi - lo > max(1, int(0.02 * hi)):
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if ok(mid) else (mid, hi)
    return {"family": family, "sd_days": sd, "n_required": hi}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--out", default=str(PROJECT_ROOT / "output" / "trackc_precision_screen"))
    args = ap.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t0, rows = time.time(), []
    for sd in SD_GRID:
        for fam in FAMILIES:
            rows.append(required_n(fam, sd, args.reps, SEED))
            print(rows[-1], flush=True)
    result = {
        "script": "src/trackc_precision_screen.py",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "rule": f"smallest n per family with {QUANTILE:.0%} quantile of GE-betting half-width <= {TARGET} days",
        "delta_per_family_contrast": DELTA, "support": SUPPORT, "reps": args.reps, "seed": SEED,
        "data": "synthetic world differences only; no policy outcomes",
        "rows": rows, "runtime_seconds": round(time.time() - t0, 1),
    }
    (out / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
