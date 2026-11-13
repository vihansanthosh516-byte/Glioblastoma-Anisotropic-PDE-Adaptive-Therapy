"""Frozen analysis for the Track C v2 evaluation run (written before any evaluation world).

For one cell: verify every planned evaluation world is complete, take each world's mean
paired difference (finite_set_minimax_mpc minus finite_set_d0_minimax_mpc) in world-index
order, and compute the GE-betting interval per response family (delta 0.025 each,
Bonferroni over 2 families = 95% familywise; support [-790, 790] days).

Rules (config outcomes_and_inference; paper/trackc_amendment_unblock_2026-11-12.md):
- positive_rule: both lower bounds > 0 AND both half-widths <= 14 days.
- An incomplete run, an empty interval (None), or a value outside the support stops
  inference with status "inconclusive"; nothing is dropped or replaced.
- Order: worlds enter the interval in schedule index order, never sorted or filtered.
- Only the primary cell can be "positive"; sensitivity cells are descriptive.
Development-pilot records are refused. Simulation only.

  python -I src/trackc_v2_analysis.py --cell primary --patients M --n-worlds N
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
from trackc_ge_betting import ge_betting_interval  # noqa: E402

PRIMARY_CONTRAST = "minimax"
SUPPORT = (-790.0, 790.0)
DELTA = 0.025
TARGET = 14.0


def family_differences(run_dir: Path, jobs, contrast=PRIMARY_CONTRAST) -> dict:
    """World means per family, in schedule order. Raises if any record is not usable."""
    out = {f: [] for f in rd.FAMILIES}
    for job in jobs:  # jobs are in schedule order
        rec = rd.read_checkpoint(rd.checkpoint_path(run_dir, job["world"]["world_id"]))
        if rec is None:
            raise ValueError(f"missing or corrupt world {job['world']['world_id']}")
        if rec["phase"] != "evaluation":
            raise ValueError("only evaluation-phase records may be analysed")
        out[job["family"]].append(float(rec["contrasts"][contrast]["mean"]))
    return out


def analyse(diffs: dict, *, primary: bool) -> dict:
    fam_rows, problems = {}, []
    for fam in rd.FAMILIES:
        x = np.asarray(diffs.get(fam, []), dtype=float)
        if x.size < 2 or not np.isfinite(x).all() or (x < SUPPORT[0]).any() or (x > SUPPORT[1]).any():
            problems.append(f"{fam}: fewer than 2 finite in-support world differences")
            fam_rows[fam] = None
            continue
        iv = ge_betting_interval(x, support=SUPPORT, delta=DELTA)
        if iv is None:
            problems.append(f"{fam}: empty GE-betting interval")
            fam_rows[fam] = None
            continue
        half = (iv[1] - iv[0]) / 2.0
        fam_rows[fam] = {"n_worlds": int(x.size), "mean_days": float(x.mean()), "lower": iv[0], "upper": iv[1],
                         "half_width": half, "meets_precision": half <= TARGET,
                         "lower_above_zero": iv[0] > 0, "upper_below_zero": iv[1] < 0}
    if problems:
        status = "inconclusive"
    else:
        rows = fam_rows.values()
        precise = all(r["meets_precision"] for r in rows)
        if precise and all(r["lower_above_zero"] for r in rows):
            status = "positive" if primary else "descriptive_both_above_zero"
        elif not precise:
            status = "inconclusive_precision_not_met"
        elif all(r["upper_below_zero"] for r in rows):
            status = "negative_both_below_zero"
        else:
            status = "not_positive"
    return {"status": status, "families": fam_rows, "problems": problems,
            "rule": "positive iff both lower bounds > 0 and both half-widths <= 14 days (primary cell only)"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", choices=rd.CELLS, default="primary")
    ap.add_argument("--patients", type=int, required=True)
    ap.add_argument("--n-worlds", type=int, required=True)
    ap.add_argument("--out", default=str(rd.DEFAULT_OUT))
    args = ap.parse_args(argv)
    jobs = rd.plan_jobs("evaluation", args.cell, args.patients, args.n_worlds, "real")
    run_dir = rd.run_dir_for(Path(args.out), "evaluation", args.cell, args.patients)
    report = rd.verify(run_dir, jobs)
    if report["status"] != "complete":
        result = {"status": "inconclusive", "problems": ["run incomplete"], "verification": report}
    else:
        result = analyse(family_differences(run_dir, jobs), primary=args.cell == "primary")
        result["verification_status"] = report["status"]
    result.update({"cell": args.cell, "patients_per_cohort": args.patients, "n_worlds_per_family": args.n_worlds,
                   "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   "label": "simulation only; not a clinical result"})
    (run_dir / "analysis.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(result["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
