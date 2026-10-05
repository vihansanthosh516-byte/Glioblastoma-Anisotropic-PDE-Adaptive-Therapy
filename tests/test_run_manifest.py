"""Run manifest writer records what reproduction needs (analysis_plan_v1.md A1.17)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.run_manifest import sha256_file, write_run_manifest  # noqa: E402


def test_manifest_fields_and_hashes(tmp_path):
    f = tmp_path / "in.txt"
    f.write_text("abc")
    out = tmp_path / "m.json"
    m = write_run_manifest("exp1", out, script="s.py", seed=7, config={"k": 1}, inputs=[f, tmp_path / "nope.txt"],
                           dataset="d", patient_split="s", primary_endpoint="e")
    disk = json.loads(out.read_text())
    assert disk == json.loads(json.dumps(m))
    for key in ["experiment_id", "seed", "config", "code_commit", "python", "packages", "inputs", "date_utc",
                "tracked_files_modified", "dataset", "patient_split", "primary_endpoint"]:
        assert key in disk
    assert len(disk["code_commit"]) == 40
    assert disk["inputs"][str(f)]["sha256"] == sha256_file(f)
    assert disk["inputs"][str(tmp_path / "nope.txt")] == {"missing": True}
    assert disk["packages"]["numpy"]


def test_hash_changes_with_content(tmp_path):
    f = tmp_path / "a"
    f.write_text("1")
    h1 = sha256_file(f)
    f.write_text("2")
    assert sha256_file(f) != h1
