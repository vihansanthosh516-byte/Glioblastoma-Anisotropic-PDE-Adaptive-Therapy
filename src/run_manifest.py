"""Run manifest writer (analysis_plan_v1.md A1.17; masterplan §104, §106).

Every primary result calls write_run_manifest() and saves the sidecar next to its output JSON.
It records what is needed to reproduce the run: seed, commit, dirty state, package versions,
hardware, config, and a SHA-256 of each input file. It never modifies the inputs.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from importlib import metadata
from pathlib import Path as _Path

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
TRACKED_PACKAGES = ("numpy", "scipy", "pandas", "scikit-learn", "torch", "torch-geometric",
                    "scanpy", "anndata", "nibabel", "lifelines", "matplotlib")


def sha256_file(path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def _git(*args: str) -> str:
    try:
        return subprocess.run(["git", *args], cwd=PROJECT_ROOT, capture_output=True, text=True,
                              timeout=30).stdout.strip()
    except Exception:
        return ""


def package_versions() -> dict:
    out = {}
    for p in TRACKED_PACKAGES:
        try:
            out[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            out[p] = None
    return out


def write_run_manifest(experiment_id: str, out_path, *, script: str, seed, config: dict | None = None,
                       inputs: list | None = None, primary_endpoint: str | None = None,
                       dataset: str | None = None, patient_split: str | None = None,
                       max_hash_bytes: int = 2_000_000_000) -> dict:
    """Write <out_path> (a .json) and return the manifest dict.

    inputs: paths of input files; each is hashed (files over max_hash_bytes record size only).
    """
    hashes = {}
    for p in inputs or []:
        p = _Path(p)
        if not p.exists():
            hashes[str(p)] = {"missing": True}
            continue
        size = p.stat().st_size
        hashes[str(p.relative_to(PROJECT_ROOT) if PROJECT_ROOT in p.parents else p)] = {
            "bytes": size, "sha256": sha256_file(p) if size <= max_hash_bytes else None}
    dirty = _git("status", "--porcelain", "--untracked-files=no")
    m = {
        "experiment_id": experiment_id,
        "script": script,
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "dataset": dataset,
        "patient_split": patient_split,
        "primary_endpoint": primary_endpoint,
        "seed": seed,
        "config": config or {},
        # inside the Docker image there is no .git; the build passes the commit in as GIT_COMMIT
        "code_commit": _git("rev-parse", "HEAD") or os.environ.get("GIT_COMMIT", ""),
        "tracked_files_modified": bool(dirty),
        "modified_files": dirty.splitlines()[:20],
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "processor": platform.processor(),
        "packages": package_versions(),
        "inputs": hashes,
    }
    out_path = _Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(m, indent=2))
    return m
