from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

from src.track_b_run_provenance import build_run_provenance


def _git(repo: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def test_run_provenance_hashes_named_inputs_and_tracks_revision_state(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "config", "user.email", "test@example.invalid")
    _git(repo, "config", "user.name", "Test")

    source = repo / "src" / "runner.py"
    manifest = repo / "data" / "manifest.csv"
    source.parent.mkdir()
    manifest.parent.mkdir()
    source.write_text("VALUE = 1\n", encoding="utf-8")
    manifest.write_text("field\nsynthetic\n", encoding="utf-8")
    _git(repo, "add", "src/runner.py", "data/manifest.csv")
    _git(repo, "commit", "-q", "-m", "synthetic fixture")
    revision = _git(repo, "rev-parse", "HEAD")

    clean = build_run_provenance(
        repo,
        manifests={"patient_manifest": manifest},
        source_files={"runner": source},
        distributions={},
    )

    assert clean["source_revision"] == revision
    assert clean["working_tree_dirty"] is False
    assert clean["manifest_sha256"]["patient_manifest"] == hashlib.sha256(
        manifest.read_bytes()
    ).hexdigest()
    assert clean["source_file_sha256"]["runner"] == hashlib.sha256(
        source.read_bytes()
    ).hexdigest()
    assert clean["python_version"]
    assert clean["package_versions"] == {}

    source.write_text("VALUE = 2\n", encoding="utf-8")
    dirty = build_run_provenance(
        repo,
        manifests={"patient_manifest": manifest},
        source_files={"runner": source},
        distributions={},
    )
    assert dirty["working_tree_dirty"] is True
    assert dirty["source_file_sha256"]["runner"] != clean["source_file_sha256"]["runner"]
