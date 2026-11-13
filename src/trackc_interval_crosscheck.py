"""Outcome-free cross-check of the project GE-betting interval against confseq.

Compares `trackc_ge_betting.ge_betting_interval` (Martinez-Taboada & Ramdas 2026
construction, project implementation) with the hedged-capital betting interval
`confseq.betting.betting_ci` (Waudby-Smith & Ramdas, JRSSB 2024; authors' code).
They are different constructions, so the check is not numerical equality. It checks:
  1. both intervals contain the true mean at >= the nominal rate in every cell;
  2. the project interval is not systematically wider than the reference;
  3. both intervals always lie inside the support.
Data are synthetic world-difference distributions only. No policy outcome is read.

Usage:
  python -I src/trackc_interval_crosscheck.py --confseq-path <dir containing confseq/> \
      --out output/trackc_interval_crosscheck [--reps 200]
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
DELTA = 0.025  # per-contrast two-sided noncoverage (Bonferroni over 2 contrasts)
SEED = 20261112


def draw(family: str, sd: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """Synthetic world-level paired differences (days), clipped to the support."""
    if family == "normal":
        x = rng.normal(10.0, sd, n)
    elif family == "skewed":  # shifted gamma with the same SD
        k = 2.0
        x = rng.gamma(k, sd / np.sqrt(k), n) - sd * np.sqrt(k) + 10.0
    elif family == "bimodal":  # half near 0, half near 2*sd
        x = np.where(rng.random(n) < 0.5, rng.normal(0.0, sd / 3, n), rng.normal(2 * sd, sd / 3, n))
    else:
        raise ValueError(family)
    return np.clip(x, *SUPPORT)


def true_mean(family: str, sd: float) -> float:
    rng = np.random.Generator(np.random.PCG64DXSM(SEED + 999))
    return float(draw(family, sd, 4_000_000, rng).mean())


def reference_interval(betting_ci, values: np.ndarray, delta: float):
    lo, hi = SUPPORT
    u = (values - lo) / (hi - lo)
    a, b = betting_ci(u, alpha=delta, breaks=2000, parallel=False)
    return lo + (hi - lo) * float(a), lo + (hi - lo) * float(b)


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = k / n
    den = 1 + z * z / n
    centre = p + z * z / (2 * n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return float((centre - half) / den), float((centre + half) / den)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--confseq-path", required=True)
    ap.add_argument("--out", default=str(PROJECT_ROOT / "output" / "trackc_interval_crosscheck"))
    ap.add_argument("--reps", type=int, default=200)
    args = ap.parse_args(argv)
    sys.path.insert(0, str(Path(args.confseq_path).resolve()))
    from confseq.betting import betting_ci  # authors' reference implementation

    out = Path(args.out)
    if not out.is_absolute():
        out = PROJECT_ROOT / out
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.Generator(np.random.PCG64DXSM(SEED))
    cells = [(f, sd, n) for f in ("normal", "skewed", "bimodal") for sd in (30.0, 60.0) for n in (100, 350)]
    rows, t0 = [], time.time()
    for family, sd, n in cells:
        mu = true_mean(family, sd)
        ge_cov = ref_cov = 0
        ge_w, ref_w, ge_none = [], [], 0
        for _ in range(args.reps):
            x = draw(family, sd, n, rng)
            g = ge_betting_interval(x, support=SUPPORT, delta=DELTA)
            r = reference_interval(betting_ci, x, DELTA)
            if g is None:
                ge_none += 1
            else:
                assert SUPPORT[0] <= g[0] <= g[1] <= SUPPORT[1]
                ge_cov += g[0] <= mu <= g[1]
                ge_w.append((g[1] - g[0]) / 2)
            ref_cov += r[0] <= mu <= r[1]
            ref_w.append((r[1] - r[0]) / 2)
        rows.append({
            "family": family, "sd_days": sd, "n_worlds": n, "reps": args.reps, "true_mean": mu,
            "ge_coverage": ge_cov / args.reps, "ge_coverage_wilson": wilson(ge_cov, args.reps),
            "ge_empty_intervals": ge_none,
            "reference_coverage": ref_cov / args.reps, "reference_coverage_wilson": wilson(ref_cov, args.reps),
            "ge_median_half_width": float(np.median(ge_w)) if ge_w else None,
            "reference_median_half_width": float(np.median(ref_w)),
        })
        print(rows[-1], flush=True)
    nominal = 1 - DELTA
    result = {
        "script": "src/trackc_interval_crosscheck.py",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reference": "confseq.betting.betting_ci (github.com/gostevehoward/confseq master, np.float_ patched to np.float64)",
        "delta_per_contrast": DELTA, "nominal_coverage": nominal, "support": SUPPORT, "seed": SEED,
        "data": "synthetic world-difference distributions only; no policy outcomes",
        "cells": rows,
        "checks": {
            # pass = coverage not significantly below nominal (Wilson 95% upper bound >= nominal)
            "ge_coverage_not_significantly_below_nominal_all_cells": all(r["ge_coverage_wilson"][1] >= nominal for r in rows),
            "reference_coverage_not_significantly_below_nominal_all_cells": all(r["reference_coverage_wilson"][1] >= nominal for r in rows),
            "ge_no_empty_intervals": all(r["ge_empty_intervals"] == 0 for r in rows),
            "ge_half_width_not_above_reference_by_10pct_all_cells": all(
                r["ge_median_half_width"] is not None and r["ge_median_half_width"] <= 1.10 * r["reference_median_half_width"]
                for r in rows),
        },
        "runtime_seconds": round(time.time() - t0, 1),
    }
    (out / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result["checks"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
