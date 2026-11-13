"""Track C integrated v2 research-scale driver (C-3, C-4, C-11, C-12, C-16).

Runs every planned world for one phase and one analysis cell, in parallel, with one
checkpoint file per world. A stopped or crashed run resumes and skips finished worlds.
At the end it verifies that every planned world has a complete, hash-checked record.

Phases
  dev_pilot   development worlds from the reserved trackc_v2_dev_pilot seed block. Used
              only to estimate variance components for the cohort-size and world-count
              choice (decision D3). Never evaluation data. Needs --confirm-development.
  evaluation  confirmatory worlds from the trackc_v2_evaluation block. Refused while any
              config gate is open (trackc_integrated_driver.check_evaluation_gates).

Cells (all settings read from configs/trackABC_integrated_v2.json at run time)
  primary, scan28, scan56, scan112, miss10, penalty0, penalty10

Inference uses the GE-betting interval only; the rejected max-T helper is refused.
Simulation only. No patient data is read.

Examples
  python -I src/trackc_research_driver.py --phase dev_pilot --cell primary --patients 12 \
      --n-worlds 20 --workers 18 --max-hours 5 --confirm-development
  python -I src/trackc_research_driver.py --phase dev_pilot --cell primary --patients 12 \
      --n-worlds 20 --verify-only
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import platform
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path

# One BLAS/OpenMP thread per worker process: the driver parallelizes over worlds, so
# multi-threaded math libraries would oversubscribe the cores (set before numpy loads).
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC = PROJECT_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import trackc_integrated_driver as integrated  # noqa: E402
import trackc_seed_registry as registry  # noqa: E402
from trackc_integrated_signal_model import F0_REFERENCE_BOUNDS, SignalModelParameters  # noqa: E402
from trackc_research_grid import research_model_grid  # noqa: E402
from trackc_world_schedule import POLICY_ROLES, SCENARIO_ROLES  # noqa: E402

CONFIG_PATH = PROJECT_ROOT / "configs" / "trackABC_integrated_v2.json"
DEFAULT_OUT = PROJECT_ROOT / "output" / "trackc_v2_research"
FAMILIES = ("single", "mixture")
CELLS = ("primary", "scan28", "scan56", "scan112", "miss10", "penalty0", "penalty10")
PHASES = {
    "dev_pilot": ("trackc_v2_dev_pilot", registry.V2_DEV_PILOT_SCENARIO_BLOCK, registry.V2_DEV_PILOT_POLICY_BLOCK),
    "evaluation": ("trackc_v2_evaluation", registry.V2_EVALUATION_SCENARIO_BLOCK, registry.V2_EVALUATION_POLICY_BLOCK),
}
INTERVAL_METHOD = "ge_betting"
REFUSED_INTERVAL_METHODS = ("simultaneous_max_t_world_intervals", "max_t")
# Signal coefficients used to generate the simulated signals (primary anchor). The world's
# q comes from the world draw (prepare_synthetic_world); a, b, d from the config anchor.
CHECKPOINT_VERSION = 1


def load_config(path=CONFIG_PATH) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _load_benchmark_module():
    name = "trackc_policy_benchmark_211"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SRC / "211_trackc_policy_benchmark.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@dataclass(frozen=True)
class CellSpec:
    name: str
    scan_every: int
    forecast_cycles: int
    course_penalty: float
    missingness: float


def cell_spec(cell: str, config=None) -> CellSpec:
    """Build one analysis cell from the config. Unknown cells are refused."""
    c = load_config() if config is None else config
    clock, ctrl, out = c["simulation_clock_and_treatment"], c["controller_comparison"], c["outcomes_and_inference"]
    cycle = int(clock["cycle_days"])
    primary_scan = int(clock["scan_interval_primary_days"]) // cycle
    forecast = int(ctrl["forecast_cycles"])
    long_forecast = int(ctrl["forecast_cycles_by_scan_interval_days"]["112"])
    penalty = _primary_penalty(out)
    miss0 = float(clock["visible_scan_missingness_primary"])
    if cell == "primary":
        return CellSpec(cell, primary_scan, forecast, penalty, miss0)
    if cell.startswith("scan"):
        days = int(cell[4:])
        if days not in clock["scan_interval_sensitivities_days"]:
            raise ValueError(f"scan interval {days} is not a declared sensitivity")
        every = days // cycle
        cycles = long_forecast if days == 112 else forecast
        return CellSpec(cell, every, cycles, penalty, miss0)
    if cell == "miss10":
        return CellSpec(cell, primary_scan, forecast, penalty, float(clock["visible_scan_missingness_sensitivity"]))
    if cell.startswith("penalty"):
        value = float(cell[7:])
        if value not in [float(v) for v in out["course_penalty_sensitivities_days"]]:
            raise ValueError(f"course penalty {value} is not a declared sensitivity")
        return CellSpec(cell, primary_scan, forecast, value, miss0)
    raise ValueError(f"unknown cell {cell!r}; choose from {CELLS}")


def _primary_penalty(outcomes) -> float:
    text = outcomes["utility"]
    marker = "minus "
    value = text.split(marker, 1)[1].split(" ", 1)[0]
    return float(value)


def make_settings(spec: CellSpec, config=None):
    """PolicyBenchmarkSettings for one cell; stop rule passed explicitly (C-3)."""
    c = load_config() if config is None else config
    bench = _load_benchmark_module()
    ctrl = c["controller_comparison"]
    return bench.PolicyBenchmarkSettings(
        scan_every=spec.scan_every,
        stop_endpoint=c["simulation_clock_and_treatment"]["stop_endpoint_code"],
        mpc_assumed_future_scan_miss_probability=spec.missingness,
        tau_days=float(c["outcomes_and_inference"]["progression_horizon_days"]),
        course_penalty=spec.course_penalty,
        mpc_particles=int(ctrl["particle_filters"]),
        prediction_particles=int(ctrl["common_predictive_draws_per_model"]),
        forecast_cycles=spec.forecast_cycles,
        signal_f0_bounds=F0_REFERENCE_BOUNDS,  # config signal_model.resistance_score_reference (parity-tested)
    )


def build_schedule(phase: str, n_worlds: int, families=FAMILIES) -> list[dict]:
    """Deterministic per-family world schedule inside the phase's reserved blocks.

    World slot g = i * len(FAMILIES) + family_index, so world i of each family keeps the
    same seeds whatever n_worlds is (worlds are nested as n grows). Every cell uses the
    same worlds; cells are separate analyses.
    """
    if phase not in PHASES:
        raise ValueError(f"phase must be one of {tuple(PHASES)}")
    if isinstance(n_worlds, bool) or not isinstance(n_worlds, int) or n_worlds < 1:
        raise ValueError("n_worlds must be a positive integer")
    purpose, scen_block, pol_block = PHASES[phase]
    ns, npol = len(SCENARIO_ROLES), len(POLICY_ROLES)
    worlds = []
    for i in range(n_worlds):
        for fi, family in enumerate(FAMILIES):
            if family not in families:
                continue
            g = i * len(FAMILIES) + fi
            w = {"world_id": f"{phase}-{family}-{i:05d}"}
            w.update({r: scen_block[0] + g * ns + j for j, r in enumerate(SCENARIO_ROLES)})
            w.update({r: pol_block[0] + g * npol + j for j, r in enumerate(POLICY_ROLES)})
            if w[SCENARIO_ROLES[-1]] > scen_block[1] or w[POLICY_ROLES[-1]] > pol_block[1]:
                raise ValueError("n_worlds exceeds the reserved seed block")
            worlds.append((family, w))
    return worlds


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, float) and not np.isfinite(obj):
        return repr(obj)
    if hasattr(obj, "__dict__") and not isinstance(obj, type):
        return _jsonable(vars(obj))
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return repr(obj)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def code_fingerprint() -> dict:
    files = ["trackc_research_driver.py", "211_trackc_policy_benchmark.py", "208_trackc_robust_belief_mpc.py",
             "trackc_world_runner.py", "trackc_integrated_sampler.py", "trackc_integrated_signal_model.py",
             "trackc_research_grid.py", "trackc_sim.py", "trackc_filter_priors.py"]
    return {f: _sha((SRC / f).read_bytes()) for f in files if (SRC / f).exists()} | {
        "config": _sha(CONFIG_PATH.read_bytes())}


def checkpoint_path(run_dir: Path, world_id: str) -> Path:
    return run_dir / "worlds" / f"{world_id}.json.gz"


def write_checkpoint(path: Path, record: dict) -> None:
    body = json.dumps(record, sort_keys=True, separators=(",", ":")).encode("utf-8")
    wrapped = json.dumps({"version": CHECKPOINT_VERSION, "sha256": _sha(body), "record": record},
                         sort_keys=True, separators=(",", ":")).encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp{os.getpid()}")
    with gzip.open(tmp, "wb") as fh:
        fh.write(wrapped)
    os.replace(tmp, path)  # atomic: a crash never leaves a half-written checkpoint


def read_checkpoint(path: Path):
    """Return the record if the file is complete and its hash matches, else None."""
    try:
        with gzip.open(path, "rb") as fh:
            doc = json.loads(fh.read().decode("utf-8"))
        body = json.dumps(doc["record"], sort_keys=True, separators=(",", ":")).encode("utf-8")
        if doc.get("version") != CHECKPOINT_VERSION or doc.get("sha256") != _sha(body):
            return None
        return doc["record"]
    except (OSError, ValueError, KeyError, EOFError):
        return None


def summarize_result(result: dict) -> dict:
    """Small per-world summary used by later analysis; the full result is kept too."""
    contrasts = result.get("paired_contrasts", {})
    summary = {}
    for key, contrast in contrasts.items():
        diff = np.asarray(contrast["robust_utility_difference"], dtype=float)
        summary[key] = {"mean": float(diff.mean()), "per_patient": diff.tolist()}
    return summary


def run_one_world(job: dict) -> dict:
    """Worker: prepare one world, run the benchmark, return a checkpoint record."""
    import trackc_world_runner as runner

    config = load_config()
    spec = CellSpec(**job["cell"])
    settings = make_settings(spec, config) if job["benchmark"] == "real" else job["stub_settings"]
    anchor = config["signal_model"]["primary_anchor"]
    signal = SignalModelParameters(q=0.5, a=float(anchor[0]), b=float(anchor[1]), d=float(anchor[2]))
    scan_loss = {"miss_p": spec.missingness} if spec.missingness > 0 else None
    t0 = time.time()
    prepared = runner.prepare_synthetic_world(
        job["world"], family=job["family"], patients_per_cohort=job["patients"],
        signal_parameters=signal, model_parameters=research_model_grid(config), settings=settings,
        scan_loss=scan_loss, seed_purpose=job["seed_purpose"],
    )
    integrated.check_arm_plan(prepared["benchmark_kwargs"])
    benchmark = _load_benchmark_module().run_shared_policy_world if job["benchmark"] == "real" else _stub_benchmark
    result = runner.execute_prepared_world(prepared, benchmark)
    full = _jsonable(result)
    return {
        "world_id": job["world"]["world_id"], "family": job["family"], "phase": job["phase"],
        "cell": job["cell"], "patients_per_cohort": job["patients"], "benchmark": job["benchmark"],
        "seeds": job["world"], "world_parameters_sha256": prepared["world_parameters_sha256"],
        "contrasts": summarize_result(result), "elapsed_seconds": round(time.time() - t0, 3),
        "host": platform.node(), "code": job["code"], "result": full,
    }


def _stub_benchmark(**kwargs):
    """Deterministic stand-in for tests: no simulation, no policy scoring."""
    n = len(kwargs["evaluation_patient_ids"])
    seed = int(kwargs["evaluation_seed"][0])
    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    diff = rng.normal(0.0, 1.0, n)
    if os.environ.get("TRACKC_STUB_FAIL") == kwargs["world_id"]:
        raise RuntimeError("stub failure requested")
    return {"world_id": kwargs["world_id"], "stub": True,
            "paired_contrasts": {"minimax": {"robust_utility_difference": diff}}}


def plan_jobs(phase, cell, patients, n_worlds, benchmark, families=FAMILIES, stub_settings=None):
    purpose = PHASES[phase][0]
    spec = cell_spec(cell)
    code = code_fingerprint()
    return [{"phase": phase, "family": fam, "world": w, "patients": patients, "seed_purpose": purpose,
             "cell": spec.__dict__, "benchmark": benchmark, "stub_settings": stub_settings, "code": code}
            for fam, w in build_schedule(phase, n_worlds, families)]


def run_dir_for(out: Path, phase: str, cell: str, patients: int) -> Path:
    return Path(out) / phase / cell / f"m{patients:02d}"


def verify(run_dir: Path, jobs) -> dict:
    missing, bad, hashes = [], [], {}
    for job in jobs:
        wid = job["world"]["world_id"]
        path = checkpoint_path(run_dir, wid)
        if not path.exists():
            missing.append(wid)
            continue
        rec = read_checkpoint(path)
        if rec is None or rec["seeds"] != job["world"] or rec["cell"] != job["cell"] \
                or rec["patients_per_cohort"] != job["patients"] or rec["code"] != job["code"]:
            bad.append(wid)
            continue
        hashes[wid] = _sha(path.read_bytes())
    return {"planned": len(jobs), "complete": len(hashes), "missing": missing, "corrupt_or_mismatched": bad,
            "status": "complete" if len(hashes) == len(jobs) else "incomplete", "checkpoint_sha256": hashes}


def execute(jobs, run_dir: Path, *, workers: int, max_hours: float | None, log=print) -> dict:
    """Run jobs not yet checkpointed. Stops submitting new worlds after max_hours."""
    run_dir.mkdir(parents=True, exist_ok=True)
    todo = []
    for j in jobs:
        rec = read_checkpoint(checkpoint_path(run_dir, j["world"]["world_id"]))
        if rec is None:
            todo.append(j)
        elif rec["code"] != j["code"]:
            raise RuntimeError(f"{j['world']['world_id']} was run with different code or config; "
                               "refusing to mix results (use a new --out folder)")
    log(f"{len(jobs) - len(todo)} of {len(jobs)} worlds already done; {len(todo)} to run")
    deadline = None if max_hours is None else time.time() + max_hours * 3600.0
    failures, done = [], 0
    if workers <= 1:
        for job in todo:
            if deadline and time.time() >= deadline:
                log("time limit reached; stopping cleanly")
                break
            try:
                write_checkpoint(checkpoint_path(run_dir, job["world"]["world_id"]), run_one_world(job))
                done += 1
            except Exception as exc:  # a failed world is recorded, never silently skipped
                failures.append({"world_id": job["world"]["world_id"], "error": repr(exc)})
                log(f"FAILED {job['world']['world_id']}: {exc!r}")
    else:
        queue = list(todo)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            running = {}
            while queue or running:
                while queue and len(running) < workers and not (deadline and time.time() >= deadline):
                    job = queue.pop(0)
                    running[pool.submit(run_one_world, job)] = job
                if not running:
                    log("time limit reached; stopping cleanly")
                    break
                finished, _ = wait(running, return_when=FIRST_COMPLETED)
                for fut in finished:
                    job = running.pop(fut)
                    wid = job["world"]["world_id"]
                    try:
                        write_checkpoint(checkpoint_path(run_dir, wid), fut.result())
                        done += 1
                        log(f"done {wid}")
                    except Exception as exc:
                        failures.append({"world_id": wid, "error": repr(exc)})
                        log(f"FAILED {wid}: {exc!r}")
    return {"ran_this_session": done, "failures": failures}


def write_seed_table(phase: str, n_worlds: int, path: Path) -> str:
    """Materialize the derived seed table (C-9). Worlds are nested, so a table for a large
    n contains every smaller run. Returns the SHA-256 of the written file."""
    rows = [{"family": fam, **w} for fam, w in build_schedule(phase, n_worlds)]
    seeds = [v for r in rows for k, v in r.items() if k not in ("family", "world_id")]
    if len(seeds) != len(set(seeds)):
        raise ValueError("seed reused inside the table")
    purpose = PHASES[phase][0]
    from trackc_world_schedule import validate_development_schedule
    for fam in FAMILIES:
        validate_development_schedule([{k: v for k, v in r.items() if k != "family"} for r in rows if r["family"] == fam],
                                      purpose=purpose)
    body = json.dumps({"phase": phase, "purpose": purpose, "n_worlds_per_family": n_worlds,
                       "derivation": "slot g = i*2 + family_index; scenario seed = block_start + 6g + role; "
                                     "policy seed = block_start + 5g + role",
                       "rows": rows}, indent=1).encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(body)
    return _sha(body)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phase", choices=tuple(PHASES), required=True)
    ap.add_argument("--cell", choices=CELLS, default="primary")
    ap.add_argument("--patients", type=int, required=True, help="virtual patients per cohort (decision D3)")
    ap.add_argument("--n-worlds", type=int, required=True, help="worlds per response family")
    ap.add_argument("--families", default="single,mixture")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--max-hours", type=float, default=None, help="stop starting new worlds after this")
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument("--confirm-development", action="store_true",
                    help="required for dev_pilot: results are development data, never evaluation data")
    ap.add_argument("--interval-method", default=INTERVAL_METHOD)
    args = ap.parse_args(argv)

    if args.interval_method in REFUSED_INTERVAL_METHODS or args.interval_method != INTERVAL_METHOD:
        raise SystemExit("refused: only the GE-betting interval is allowed; the max-T method failed its screen")
    if args.phase == "evaluation":
        integrated.check_evaluation_gates(integrated.load_protocol())  # raises while any gate is open
    elif not args.verify_only and not args.confirm_development:
        raise SystemExit("dev_pilot needs --confirm-development (development data only)")
    families = tuple(f for f in args.families.split(",") if f)
    if not set(families) <= set(FAMILIES):
        raise SystemExit(f"families must be from {FAMILIES}")

    jobs = plan_jobs(args.phase, args.cell, args.patients, args.n_worlds, "real", families)
    run_dir = run_dir_for(Path(args.out), args.phase, args.cell, args.patients)
    session = None
    if not args.verify_only:
        session = execute(jobs, run_dir, workers=args.workers, max_hours=args.max_hours)
    report = verify(run_dir, jobs)
    manifest = {
        "phase": args.phase, "cell": cell_spec(args.cell).__dict__, "patients_per_cohort": args.patients,
        "n_worlds_per_family": args.n_worlds, "families": list(families), "interval_method": INTERVAL_METHOD,
        "model_count": len(research_model_grid()), "code": code_fingerprint(), "host": platform.node(),
        "python": sys.version, "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "session": session, "verification": report,
        "label": "development data; never evaluation data" if args.phase == "dev_pilot" else "evaluation",
    }
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"{report['status']}: {report['complete']} of {report['planned']} worlds; "
          f"{len(report['missing'])} missing, {len(report['corrupt_or_mismatched'])} bad")
    return 0 if report["status"] == "complete" else 2


if __name__ == "__main__":
    raise SystemExit(main())
