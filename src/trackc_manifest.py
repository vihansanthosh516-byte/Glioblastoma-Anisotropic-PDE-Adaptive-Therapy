"""Run manifest and protocol hash for Track C runs (SIMULATION bookkeeping; runs no simulation).

Why: the 2026-10 sealed protocols hashed only the protocol text. The solver mode changed from
"grid" to "refine" at script 167 without a hash change. Here the protocol hash covers the
protocol, the solver constants (EVENT_MODE, DT, PROG_FACTOR, RANGES and the cycle constants) and
the SHA-256 of the core code files, so any of them changing changes the hash.

A manifest also records the seed blocks (checked against src/trackc_seed_registry.py), the git
commit and whether the core code files are committed. After a run that spends a reserved seed
block, move that block from "reserved" to "dev_used" or "spent" in the registry and commit.

Use:
    from trackc_manifest import build_manifest, write_manifest
    m = build_manifest("stageA_tail", protocol, purpose="stage2_dev",
                       scenario=(20261800, 20261899), noise=(1800000, 1800999),
                       script="src/<your_script>.py")
    write_manifest(OUT_DIR, m)           # refuses to overwrite an existing manifest.json
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import trackc_seed_registry as reg
import trackc_sim as tc

ROOT = Path(__file__).resolve().parent.parent
CORE_CODE = ("src/trackc_sim.py", "src/trackc_oracle.py", "src/trackc_reference.py",
             "src/trackc_seed_registry.py", "src/trackc_manifest.py")
LABEL = "SIMULATION ONLY; manifest of a Track C run; no clinical meaning"


def solver_constants():
    """Module constants that change simulated outcomes. Read at call time, so a patched value shows."""
    return {
        "EVENT_MODE": tc.EVENT_MODE, "DT": tc.DT, "PROG_FACTOR": tc.PROG_FACTOR,
        "CYCLE_DAYS": tc.CYCLE_DAYS, "COURSE_DAYS": tc.COURSE_DAYS, "N_CYCLES": tc.N_CYCLES,
        "HORIZON": tc.HORIZON, "BUDGET": tc.BUDGET, "SIGMAS": list(tc.SIGMAS),
        "RANGES": {k: list(v) for k, v in tc.RANGES.items()},
        "RANGES_WIDE_M": {k: list(v) for k, v in tc.RANGES_WIDE_M.items()},
    }


def code_hashes(paths=CORE_CODE, root=ROOT):
    """SHA-256 of each code file, or None if the file is missing."""
    out = {}
    for rel in paths:
        p = Path(root) / rel
        out[rel] = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
    return out


def protocol_hash(protocol, paths=CORE_CODE, root=ROOT):
    """SHA-256 over the protocol, the solver constants and the code hashes."""
    body = {"protocol": protocol, "solver_constants": solver_constants(),
            "code_sha256": code_hashes(paths, root)}
    return hashlib.sha256(json.dumps(body, sort_keys=True, default=str).encode("utf-8")).hexdigest()


def _git(args, root):
    try:
        res = subprocess.run(["git", *args], cwd=str(root), capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return res.stdout if res.returncode == 0 else None


def git_state(paths=CORE_CODE, root=ROOT):
    """Commit and whether the given paths have uncommitted or untracked changes.
    Only these paths are checked, because the rest of the repository is often dirty."""
    head = _git(["rev-parse", "HEAD"], root)
    status = _git(["status", "--porcelain", "--", *paths], root)
    if head is None or status is None:
        return {"commit": None, "dirty_paths": None, "clean": None}
    dirty = [line[3:].strip() for line in status.splitlines() if line.strip()]
    return {"commit": head.strip(), "dirty_paths": dirty, "clean": not dirty}


def build_manifest(name, protocol, purpose, scenario, noise, script, require_clean=True,
                   paths=CORE_CODE, root=ROOT):
    """Return the manifest dict. Raises ValueError if a seed block is not free for this purpose,
    and RuntimeError if the code files are not committed (when require_clean is True)."""
    problems = reg.check_free(scenario=scenario, noise=noise, purpose=purpose)
    if problems:
        raise ValueError("seed problems: " + "; ".join(problems))
    run_paths = tuple(paths) + ((script,) if script else ())
    state = git_state(run_paths, root)
    if require_clean and state["clean"] is not True:
        raise RuntimeError(f"code is not committed or git is unavailable: {state}")
    hashes = code_hashes(run_paths, root)
    missing = [p for p, h in hashes.items() if h is None]
    if missing:
        raise FileNotFoundError(f"missing code files: {missing}")
    return {
        "name": name, "label": LABEL, "purpose": purpose,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scenario_seeds": list(scenario), "noise_keys": list(noise) if noise is not None else [],
        "registry_sha256": reg.registry_document()["sha256"],
        "protocol": protocol, "protocol_sha256": protocol_hash(protocol, run_paths, root),
        "solver_constants": solver_constants(), "code_sha256": hashes, "git": state,
        "python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform(),
    }


def write_manifest(out_dir, manifest):
    """Write manifest.json into out_dir. Refuses to overwrite."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "manifest.json"
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.write_text(json.dumps(manifest, indent=2, default=str) + "\n", encoding="utf-8")
    return path
