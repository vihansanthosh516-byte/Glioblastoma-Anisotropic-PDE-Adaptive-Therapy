"""Choose cohort size m and world count n from the development pilot (decision D3).

Reads development-pilot checkpoints written by trackc_research_driver (phase dev_pilot),
estimates the variance of world-level paired differences as tau^2 + sigma^2 / m, and
picks the m in M_GRID that minimizes predicted total compute subject to the 14-day
precision rule from the outcome-blind screen table (src/trackc_precision_screen.py).

Rules (fixed before the pilot runs; paper/trackc_amendment_unblock_2026-11-12.md):
- sigma^2: pooled within-world variance of per-patient primary differences.
- tau^2: variance of world means minus sigma^2 / m0, floored at 0.
- Uncertainty: 2,000 bootstrap resamples of pilot worlds; use the 80th percentile of
  the predicted variance for each m (errs toward more worlds).
- n(m): screen table, worst family shape, linear interpolation in SD; round up.
- Time per world: a + b x m, fitted from the mean seconds of two pilot runs at different
  cohort sizes (--timing-dir). With one run only, time is assumed proportional to m, which
  always favours the smallest m (n(m) x m grows with m when n is proportional to variance);
  the output then says so. Optimal-cluster-size theory gives m* = sqrt(sigma^2 a / (tau^2 b)).
- Choose m with the smallest 2 x n(m) x time(m); a tie goes to the larger m.
Development data only. The output is a design choice, never an evaluation result.

  python -I src/trackc_dev_pilot_design.py --pilot-dir output/trackc_v2_research/dev_pilot/primary/m12
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_research_driver as rd  # noqa: E402

M_GRID = (2, 4, 6, 8, 12, 16)
BOOTSTRAP = 2000
UPPER_QUANTILE = 0.80
SEED = 20261114
SCREEN = PROJECT_ROOT / "output" / "trackc_precision_screen" / "result.json"


def load_pilot(pilot_dir: Path) -> list[dict]:
    recs = []
    for path in sorted((pilot_dir / "worlds").glob("*.json.gz")):
        rec = rd.read_checkpoint(path)
        if rec is None:
            raise ValueError(f"corrupt checkpoint {path}")
        if rec["phase"] != "dev_pilot":
            raise ValueError("only development-pilot records may be used")
        recs.append(rec)
    if len(recs) < 4:
        raise ValueError("need at least 4 pilot worlds")
    return recs


def variance_components(per_world: list[np.ndarray]) -> tuple[float, float, int]:
    m0 = len(per_world[0])
    if any(len(x) != m0 for x in per_world) or m0 < 2:
        raise ValueError("all pilot worlds need the same cohort size >= 2")
    sigma2 = float(np.mean([np.var(x, ddof=1) for x in per_world]))
    means = np.array([x.mean() for x in per_world])
    tau2 = max(0.0, float(np.var(means, ddof=1)) - sigma2 / m0)
    return tau2, sigma2, m0


def n_for_sd(sd: float, screen: dict) -> int | None:
    """Worst family shape; linear interpolation of n in SD; None if above the screen."""
    best = None
    for fam in {r["family"] for r in screen["rows"]}:
        rows = sorted((r["sd_days"], r["n_required"]) for r in screen["rows"] if r["family"] == fam)
        sds = [s for s, n in rows if n is not None]
        ns = [n for s, n in rows if n is not None]
        if sd > sds[-1]:
            return None
        n = int(np.ceil(np.interp(sd, sds, ns)))
        best = n if best is None else max(best, n)
    return best


def design(recs: list[dict], screen: dict, contrast="minimax", timing=None) -> dict:
    """timing: None or (a_seconds, b_seconds_per_patient) for time per world = a + b m."""
    out = {}
    rng = np.random.Generator(np.random.PCG64DXSM(SEED))
    secs = float(np.mean([r["elapsed_seconds"] for r in recs]))
    per_family = {}
    for fam in sorted({r["family"] for r in recs}):
        worlds = [np.asarray(r["contrasts"][contrast]["per_patient"], float) for r in recs if r["family"] == fam]
        tau2, sigma2, m0 = variance_components(worlds)
        boot = []
        for _ in range(BOOTSTRAP):
            idx = rng.integers(0, len(worlds), len(worlds))
            t, s, _m = variance_components([worlds[i] for i in idx])
            boot.append((t, s))
        boot = np.array(boot)
        per_family[fam] = {"n_pilot_worlds": len(worlds), "m0": m0, "tau2": tau2, "sigma2": sigma2,
                           "var_upper": {m: float(np.quantile(boot[:, 0] + boot[:, 1] / m, UPPER_QUANTILE))
                                         for m in M_GRID}}
    m0 = recs[0]["patients_per_cohort"]
    a, b = (0.0, secs / m0) if timing is None else timing
    for m in M_GRID:
        sd = max(np.sqrt(per_family[f]["var_upper"][m]) for f in per_family)
        n = n_for_sd(float(sd), screen)
        hours = None if n is None else 2 * n * (a + b * m) / 3600.0
        out[m] = {"sd_upper_days": float(sd), "n_per_family": n, "single_core_hours": hours}
    feasible = {m: v for m, v in out.items() if v["single_core_hours"] is not None}
    choice = min(feasible, key=lambda m: (feasible[m]["single_core_hours"], -m)) if feasible else None
    return {"per_family": per_family, "by_m": out, "choice_m": choice,
            "choice_n_per_family": None if choice is None else out[choice]["n_per_family"],
            "pilot_mean_seconds_per_world": secs, "time_model_seconds": {"a": a, "b": b},
            "time_model_source": "assumed proportional to m (one run)" if timing is None else "fitted from two cohort sizes",
            "label": "development design choice; not an evaluation result"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot-dir", required=True)
    ap.add_argument("--screen", default=str(SCREEN))
    ap.add_argument("--timing-dir", default=None, help="second pilot run at another cohort size, for a + b m")
    args = ap.parse_args(argv)
    pilot = Path(args.pilot_dir)
    if not pilot.is_absolute():
        pilot = PROJECT_ROOT / pilot
    screen = json.loads(Path(args.screen).read_text(encoding="utf-8"))
    recs = load_pilot(pilot)
    timing = None
    if args.timing_dir:
        other = load_pilot(Path(args.timing_dir) if Path(args.timing_dir).is_absolute() else PROJECT_ROOT / args.timing_dir)
        (m1, t1), (m2, t2) = [(r[0]["patients_per_cohort"], float(np.mean([x["elapsed_seconds"] for x in r])))
                              for r in (recs, other)]
        if m1 == m2:
            raise ValueError("timing run needs a different cohort size")
        b = (t1 - t2) / (m1 - m2)
        timing = (max(0.0, t1 - b * m1), max(b, 1e-9))
    result = design(recs, screen, timing=timing)
    result["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (pilot / "design.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("choice_m", "choice_n_per_family")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
