"""Synthetic tests for the research driver. Stub benchmark only; no policy is scored."""
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

import trackc_research_driver as rd
import trackc_research_grid as grid
import trackc_seed_registry as registry

STUB_SETTINGS = SimpleNamespace(stop_endpoint="observed_two_scan")  # picklable for worker processes


def jobs(n=2, families=rd.FAMILIES):
    return rd.plan_jobs("dev_pilot", "primary", 2, n, "stub", families, stub_settings=STUB_SETTINGS)


def test_research_grid_has_135_models_covering_every_true_q():
    models = grid.research_model_grid()
    assert len(models) == 135 == len(set(models))
    assert {0.6, 0.7, 0.8, 0.9} <= {m.q for m in models}


def test_cells_read_from_config():
    assert rd.cell_spec("primary") == rd.CellSpec("primary", 3, 6, 5.0, 0.0)
    assert rd.cell_spec("scan112") == rd.CellSpec("scan112", 4, 8, 5.0, 0.0)   # decision D2
    assert rd.cell_spec("scan28").scan_every == 1
    assert rd.cell_spec("miss10").missingness == 0.1
    assert rd.cell_spec("penalty10").course_penalty == 10.0
    with pytest.raises(ValueError):
        rd.cell_spec("scan99")


def test_every_cell_builds_valid_settings_with_explicit_stop_rule():
    for cell in rd.CELLS:
        s = rd.make_settings(rd.cell_spec(cell))
        assert s.stop_endpoint == "observed_two_scan"
        assert s.forecast_cycles % s.scan_every == 0
        assert s.signal_f0_bounds == rd.F0_REFERENCE_BOUNDS


def test_schedule_is_nested_inside_reserved_blocks_and_free():
    small, big = rd.build_schedule("dev_pilot", 2), rd.build_schedule("dev_pilot", 5)
    assert big[:len(small)] == small
    seeds = [v for _, w in big for k, v in w.items() if k != "world_id"]
    assert len(seeds) == len(set(seeds))
    lo, hi = registry.V2_DEV_PILOT_SCENARIO_BLOCK
    assert all(lo <= w["world"] <= hi for _, w in big)
    ev = rd.build_schedule("evaluation", 3)
    assert not {v for _, w in ev for k, v in w.items() if k != "world_id"} & set(seeds)


def test_end_to_end_stub_run_and_verify(tmp_path):
    js = jobs(2)
    out = rd.execute(js, tmp_path, workers=1, max_hours=None, log=lambda *_: None)
    assert out == {"ran_this_session": 4, "failures": []}
    rep = rd.verify(tmp_path, js)
    assert rep["status"] == "complete" and rep["complete"] == 4


def test_resume_skips_finished_worlds(tmp_path):
    js = jobs(2)
    rd.execute(js[:2], tmp_path, workers=1, max_hours=None, log=lambda *_: None)
    out = rd.execute(js, tmp_path, workers=1, max_hours=None, log=lambda *_: None)
    assert out["ran_this_session"] == 2


def test_verify_catches_missing_and_corrupt_world(tmp_path):
    js = jobs(2)
    rd.execute(js, tmp_path, workers=1, max_hours=None, log=lambda *_: None)
    paths = sorted((tmp_path / "worlds").glob("*.json.gz"))
    paths[0].unlink()
    paths[1].write_bytes(b"garbage")
    rep = rd.verify(tmp_path, js)
    assert rep["status"] == "incomplete"
    assert len(rep["missing"]) == 1 and len(rep["corrupt_or_mismatched"]) == 1


def test_failed_world_is_recorded_not_skipped(tmp_path, monkeypatch):
    js = jobs(1)
    monkeypatch.setenv("TRACKC_STUB_FAIL", js[0]["world"]["world_id"])
    out = rd.execute(js, tmp_path, workers=1, max_hours=None, log=lambda *_: None)
    assert len(out["failures"]) == 1
    assert rd.verify(tmp_path, js)["status"] == "incomplete"


def test_time_limit_stops_cleanly(tmp_path):
    out = rd.execute(jobs(2), tmp_path, workers=1, max_hours=0.0, log=lambda *_: None)
    assert out["ran_this_session"] == 0


def test_parallel_workers_match_serial(tmp_path):
    js = jobs(2)
    rd.execute(js, tmp_path / "a", workers=1, max_hours=None, log=lambda *_: None)
    rd.execute(js, tmp_path / "b", workers=2, max_hours=None, log=lambda *_: None)
    for j in js:
        wid = j["world"]["world_id"]
        a = rd.read_checkpoint(rd.checkpoint_path(tmp_path / "a", wid))
        b = rd.read_checkpoint(rd.checkpoint_path(tmp_path / "b", wid))
        assert a["contrasts"] == b["contrasts"]


def test_cli_refuses_evaluation_max_t_and_unconfirmed_pilot(tmp_path):
    import trackc_integrated_driver as integrated
    with pytest.raises(integrated.EvaluationBlocked):
        rd.main(["--phase", "evaluation", "--patients", "2", "--n-worlds", "1", "--out", str(tmp_path)])
    with pytest.raises(SystemExit, match="GE-betting"):
        rd.main(["--phase", "dev_pilot", "--patients", "2", "--n-worlds", "1", "--interval-method", "max_t"])
    with pytest.raises(SystemExit, match="confirm-development"):
        rd.main(["--phase", "dev_pilot", "--patients", "2", "--n-worlds", "1", "--out", str(tmp_path)])


def test_driver_does_not_import_rejected_max_t_helper():
    text = Path(rd.__file__).read_text(encoding="utf-8")
    assert "import trackc_world_summary" not in text and "from trackc_world_summary" not in text


def test_resume_refuses_checkpoints_from_changed_code(tmp_path):
    js = jobs(1)
    rd.execute(js, tmp_path, workers=1, max_hours=None, log=lambda *_: None)
    changed = [dict(j, code={**j["code"], "config": "different"}) for j in js]
    with pytest.raises(RuntimeError, match="different code"):
        rd.execute(changed, tmp_path, workers=1, max_hours=None, log=lambda *_: None)
    assert rd.verify(tmp_path, changed)["status"] == "incomplete"


def test_seed_table_is_written_and_validated(tmp_path):
    sha = rd.write_seed_table("evaluation", 50, tmp_path / "t.json")
    assert len(sha) == 64 and (tmp_path / "t.json").stat().st_size > 0


def test_worker_math_libraries_are_single_threaded():
    assert os.environ.get("OMP_NUM_THREADS") == "1"
