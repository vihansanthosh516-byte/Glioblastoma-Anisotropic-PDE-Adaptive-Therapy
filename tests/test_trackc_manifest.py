"""Track C run manifest (src/trackc_manifest.py): the protocol hash moves when the protocol, a solver
constant or a code file moves; the clean-tree check sees uncommitted changes; seeds are checked
against the registry; manifests are never overwritten. Bookkeeping only; no simulation runs."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_manifest as man  # noqa: E402
import trackc_sim as tc  # noqa: E402

PROTO = {"hypothesis": "demo", "margin_days": -7}


def _git_available():
    return shutil.which("git") is not None


def _make_repo(tmp_path, files):
    """Tiny git repo holding copies of the named code files."""
    for rel in files:
        dst = tmp_path / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(f"# {rel}\n", encoding="utf-8")
    run = lambda *a: subprocess.run(["git", *a], cwd=tmp_path, check=True, capture_output=True)
    run("init", "-q")
    run("config", "user.email", "t@example.org")
    run("config", "user.name", "t")
    run("config", "commit.gpgsign", "false")
    run("add", "-A")
    run("commit", "-q", "-m", "init")


def test_solver_constants_cover_the_things_that_changed_before():
    c = man.solver_constants()
    assert c["EVENT_MODE"] == "refine" and c["DT"] == 0.25 and c["PROG_FACTOR"] == 1.2
    assert c["BUDGET"] == 12 and c["N_CYCLES"] == 52 and c["RANGES"]["m"] == [6.0, 20.0, "log"]


def test_protocol_hash_is_stable_and_sensitive():
    base = man.protocol_hash(PROTO)
    assert man.protocol_hash(PROTO) == base and len(base) == 64
    assert man.protocol_hash(dict(PROTO, margin_days=-8)) != base


def test_protocol_hash_moves_with_solver_constants(monkeypatch):
    base = man.protocol_hash(PROTO)
    monkeypatch.setattr(tc, "EVENT_MODE", "grid")
    assert man.protocol_hash(PROTO) != base
    monkeypatch.undo()
    monkeypatch.setattr(tc, "DT", 0.1)
    assert man.protocol_hash(PROTO) != base
    monkeypatch.undo()
    assert man.protocol_hash(PROTO) == base


def test_protocol_hash_moves_with_code(tmp_path):
    files = ("src/trackc_sim.py", "src/a.py")
    _make_repo(tmp_path, files)
    h1 = man.protocol_hash(PROTO, files, tmp_path)
    (tmp_path / "src" / "a.py").write_text("# changed\n", encoding="utf-8")
    assert man.protocol_hash(PROTO, files, tmp_path) != h1


@pytest.mark.skipif(not _git_available(), reason="git not installed")
def test_git_state_detects_dirty_and_untracked(tmp_path):
    files = ("src/trackc_sim.py", "src/a.py")
    _make_repo(tmp_path, files)
    s = man.git_state(files, tmp_path)
    assert s["clean"] is True and s["dirty_paths"] == [] and len(s["commit"]) == 40
    (tmp_path / "src" / "a.py").write_text("# edited\n", encoding="utf-8")
    s = man.git_state(files, tmp_path)
    assert s["clean"] is False and "src/a.py" in s["dirty_paths"]
    # a change outside the checked paths does not count
    (tmp_path / "src" / "a.py").write_text("# src/a.py\n", encoding="utf-8")
    (tmp_path / "other.txt").write_text("x", encoding="utf-8")
    assert man.git_state(files, tmp_path)["clean"] is True


def test_git_state_without_git_repo_is_unknown(tmp_path):
    s = man.git_state(("src/x.py",), tmp_path)
    assert s == {"commit": None, "dirty_paths": None, "clean": None}


@pytest.mark.skipif(not _git_available(), reason="git not installed")
def test_build_manifest_checks_seeds_and_tree(tmp_path):
    files = ("src/trackc_sim.py", "src/run_demo.py")
    _make_repo(tmp_path, files)
    kw = dict(name="demo", protocol=PROTO, script="src/run_demo.py", paths=("src/trackc_sim.py",), root=tmp_path)
    m = man.build_manifest(purpose="stage2_dev", scenario=(20261900, 20261949), noise=(1850000, 1850499), **kw)
    assert m["git"]["clean"] is True and m["protocol_sha256"] and m["purpose"] == "stage2_dev"
    assert set(m["code_sha256"]) == {"src/trackc_sim.py", "src/run_demo.py"}
    assert m["label"].startswith("SIMULATION ONLY")
    # spent seeds are refused
    with pytest.raises(ValueError):
        man.build_manifest(purpose="stage2_dev", scenario=(20261122, 20261122), noise=(1850000, 1850001), **kw)
    # the sealed block cannot be spent under the dev purpose
    with pytest.raises(ValueError):
        man.build_manifest(purpose="stage2_dev", scenario=(20262300, 20262310), noise=(2300000, 2300010), **kw)
    # uncommitted code is refused unless asked otherwise
    (tmp_path / "src" / "run_demo.py").write_text("# edited\n", encoding="utf-8")
    with pytest.raises(RuntimeError):
        man.build_manifest(purpose="stage2_dev", scenario=(20261900, 20261949), noise=(1850000, 1850499), **kw)
    m2 = man.build_manifest(purpose="stage2_dev", scenario=(20261900, 20261949), noise=(1850000, 1850499),
                            require_clean=False, **kw)
    assert m2["git"]["clean"] is False and m2["git"]["dirty_paths"] == ["src/run_demo.py"]


def test_build_manifest_refuses_missing_code(tmp_path):
    with pytest.raises((FileNotFoundError, RuntimeError)):
        man.build_manifest("demo", PROTO, "stage2_dev", (20261900, 20261949), (1850000, 1850499),
                           script="src/nope.py", require_clean=False, paths=("src/also_nope.py",), root=tmp_path)


def test_write_manifest_never_overwrites(tmp_path):
    path = man.write_manifest(tmp_path / "run", {"name": "demo", "when": "now"})
    assert json.loads(path.read_text(encoding="utf-8"))["name"] == "demo"
    with pytest.raises(FileExistsError):
        man.write_manifest(tmp_path / "run", {"name": "again"})
