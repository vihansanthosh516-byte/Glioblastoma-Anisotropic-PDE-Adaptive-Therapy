"""Synthetic tests for the Track B pre-run seal guard. No patient data."""
import json
from dataclasses import replace

import numpy as np
import pytest

from src import fk_lumiere_result_writer as rw
from src import fk_lumiere_runner as lr
from src import fk_track_b_run_guard as g
from tests.test_fk_lumiere_protocol import _spread_cases
from tests.test_fk_lumiere_runner import TINY, FakeProvider, _folds

N_BOOT = 200


def _factory(rho_b, d_b):
    return FakeProvider()


@pytest.fixture(autouse=True)
def _clean_git(monkeypatch):
    monkeypatch.setattr(g, "_git_probe", lambda root, files: {"head": "a" * 40, "clean": True})


@pytest.fixture
def setup(tmp_path):
    cases = _spread_cases()
    folds = _folds(cases)
    seal = g.build_seal(cases=cases, fold_by_patient=folds,
                        provider_factory=_factory, recipe=TINY, n_bootstrap=N_BOOT)
    path = g.write_seal(seal, tmp_path / "seal.json")
    return cases, folds, seal, path


def _call(cases, folds, path, **kw):
    kw.setdefault("provider_factory", _factory)
    kw.setdefault("recipe", TINY)
    kw.setdefault("n_bootstrap", N_BOOT)
    kw.setdefault("allow_synthetic", True)
    if "expected_seal_sha256" not in kw:
        kw["expected_seal_sha256"] = g.file_sha256(path) if path is not None and path.is_file() else "0" * 64
    return g.run_sealed_protocol(cases, folds, seal_path=path, **kw)


def _tamper(path, fn):
    s = json.loads(path.read_text())
    fn(s)
    path.write_text(json.dumps(s))


def test_matching_seal_passes_and_runs_protocol(setup):
    cases, folds, _, path = setup
    out = _call(cases, folds, path)
    assert out["status"] == "completed (synthetic)"


def test_seal_records_all_required_fields(setup):
    _, _, seal, _ = setup
    assert set(seal["provider"]) == {"primary", "one_decade_narrower", "one_decade_wider"}
    assert seal["recipe"]["epochs"] == 1 and seal["n_bootstrap"] == N_BOOT
    assert len(seal["fold_manifest_sha256"]) == 64
    assert any(f.endswith("fk_lumiere_protocol.py") for f in seal["code_sha256"])
    assert any("amendment" in f for f in seal["protocol_sha256"])


def test_write_seal_refuses_overwrite(setup):
    _, _, seal, path = setup
    with pytest.raises(FileExistsError):
        g.write_seal(seal, path)


def test_no_seal_refused(setup, tmp_path):
    cases, folds, _, _ = setup
    with pytest.raises(g.SealMismatch, match="no seal"):
        _call(cases, folds, tmp_path / "missing.json")
    with pytest.raises(g.SealMismatch, match="no seal"):
        _call(cases, folds, None)


def test_real_cases_without_seal_refused_before_any_target_read(tmp_path):
    real, reads = _real(_spread_cases())
    with pytest.raises(g.SealMismatch):
        _call(real, _folds(real), tmp_path / "missing.json", allow_synthetic=False, gates={x: True for x in lr.LUMIERE_GATES})
    assert reads == []


def test_real_cases_with_mismatching_seal_refused_before_target_read(setup):
    _, _, _, path = setup
    real, reads = _real(_spread_cases())
    with pytest.raises(g.SealMismatch):
        _call(real, _folds(real), path, recipe=replace(TINY, epochs=2), allow_synthetic=False,
              gates={x: True for x in lr.LUMIERE_GATES})
    assert reads == []


def test_different_provider_class_refused(setup):
    cases, folds, _, path = setup

    class Other(FakeProvider):
        pass
    with pytest.raises(g.SealMismatch, match="provider"):
        _call(cases, folds, path, provider_factory=lambda r, d: Other())


def test_different_provider_settings_refused(setup):
    cases, folds, _, path = setup

    def tweaked(r, d):
        p = FakeProvider()
        p.n_draws = 7
        return p
    with pytest.raises(g.SealMismatch, match="provider"):
        _call(cases, folds, path, provider_factory=tweaked)


def test_different_recipe_refused(setup):
    cases, folds, _, path = setup
    with pytest.raises(g.SealMismatch, match="recipe"):
        _call(cases, folds, path, recipe=replace(TINY, learning_rate=0.5))


def test_different_bootstrap_refused(setup):
    cases, folds, _, path = setup
    with pytest.raises(g.SealMismatch, match="bootstrap"):
        _call(cases, folds, path, n_bootstrap=N_BOOT + 1)
    with pytest.raises(g.SealMismatch, match="bootstrap"):
        _call(cases, folds, path, n_bootstrap=None)  # default 10,000 differs from the sealed 200


def test_different_folds_refused(setup):
    cases, folds, _, path = setup
    moved = dict(folds)
    first = cases[0].patient_id
    moved[first] = (moved[first] + 1) % 5
    with pytest.raises(g.SealMismatch, match="fold manifest"):
        _call(cases, moved, path)


def test_different_patient_set_refused(setup):
    cases, folds, _, path = setup
    with pytest.raises(g.SealMismatch, match="patient set"):
        _call(cases[:-1], folds, path)


def test_changed_code_hash_refused(setup):
    cases, folds, _, path = setup
    _tamper(path, lambda s: s["code_sha256"].update({"src/fk_lumiere_runner.py": "0" * 64}))
    with pytest.raises(g.SealMismatch, match="code: file changed src/fk_lumiere_runner.py"):
        _call(cases, folds, path)


def test_added_or_removed_code_file_refused(setup):
    cases, folds, _, path = setup
    _tamper(path, lambda s: s["code_sha256"].pop("src/fk_lumiere_stats.py"))
    with pytest.raises(g.SealMismatch, match="not in seal"):
        _call(cases, folds, path)


def test_changed_protocol_hash_refused(setup):
    cases, folds, _, path = setup
    key = sorted(json.loads(path.read_text())["protocol_sha256"])[0]
    _tamper(path, lambda s: s["protocol_sha256"].update({key: "f" * 64}))
    with pytest.raises(g.SealMismatch, match="protocol: file changed"):
        _call(cases, folds, path)


def test_malformed_seal_refused(setup):
    cases, folds, _, path = setup
    path.write_text("{not json")
    with pytest.raises(g.SealMismatch, match="not valid JSON"):
        _call(cases, folds, path)
    path.write_text(json.dumps({"provider": {}}))
    with pytest.raises(g.SealMismatch, match="missing fields"):
        _call(cases, folds, path)


def test_fold_manifest_hash_is_order_independent_and_sensitive():
    a = g.fold_manifest_sha256({"a": 0, "b": 1}, ["b", "a"])
    assert a == g.fold_manifest_sha256({"b": 1, "a": 0}, ["a", "b"])
    assert a != g.fold_manifest_sha256({"a": 1, "b": 1}, ["a", "b"])


def test_default_describe_uses_real_provider_class_without_running_it():
    from src import fk_lumiere_protocol as lp
    d = g.describe_provider_factory(lp.default_provider_factory)
    assert d["primary"]["class"].endswith("LumierePDEProvider")
    assert d["primary"]["attributes"]["settings"]["rho_bounds"] == [1e-5, 1e-1]
    assert d["one_decade_wider"]["attributes"]["settings"]["D_bounds"] == [1e-6, 100.0]
    assert np.isfinite(d["primary"]["attributes"]["seed"])


ALL_GATES = {x: True for x in lr.LUMIERE_GATES}


def _real(cases):
    reads = []
    return [replace(c, synthetic=False, load_target=lambda c=c: reads.append(c.patient_id),
                    qc_record=_qc(c.qc_passed)) for c in cases], reads


def _qc(passed):
    return {"passed": bool(passed), "reasons": [] if passed else ["synthetic failure"],
            "max_rotation_deg": 0.0, "shift_mm": 0.0}


def test_real_case_without_token_raises_in_runner_and_protocol():
    from src import fk_lumiere_protocol as lp
    real, reads = _real(_spread_cases())
    with pytest.raises(PermissionError, match="run_sealed_protocol"):
        lr.run_lumiere_cv(real, _folds(real), pde_provider=FakeProvider(), gates=ALL_GATES, recipe=TINY,
                          n_bootstrap=N_BOOT)
    with pytest.raises(PermissionError, match="run_sealed_protocol"):
        lp.run_lumiere_protocol(real, _folds(real), provider_factory=_factory, gates=ALL_GATES, recipe=TINY,
                                n_bootstrap=N_BOOT)
    with pytest.raises(PermissionError, match="run_sealed_protocol"):
        lr.run_lumiere_cv(real, _folds(real), pde_provider=FakeProvider(), gates=ALL_GATES, recipe=TINY,
                          n_bootstrap=N_BOOT, _sealed_token=object())
    assert reads == []


def test_guard_rejects_synthetic_cases_in_real_entry(setup):
    cases, folds, _, path = setup
    with pytest.raises(PermissionError, match="synthetic=False"):
        _call(cases, folds, path, allow_synthetic=False)


def test_guard_rejects_false_or_missing_gates(setup):
    _, folds, _, path = setup
    real, reads = _real(_spread_cases())
    bad = dict(ALL_GATES, scoring_authorized_by_user=False)
    for gates in (bad, {}, None):
        with pytest.raises(PermissionError, match="gates not recorded"):
            _call(real, folds, path, allow_synthetic=False, gates=gates)
    assert reads == []


def test_guard_membership_in_execution_code_files():
    files = rw.execution_code_files()
    assert "src/fk_track_b_run_guard.py" in files and "src/fk_two_history_shape.py" in files
    assert "src/fk_track_b_run_guard.py" in rw.REQUIRED_CODE_FILES
    assert "src/fk_two_history_shape.py" in rw.REQUIRED_CODE_FILES


def test_seal_records_git_versions_and_manifest(setup):
    _, _, seal, _ = setup
    assert seal["git"] == {"head": "a" * 40, "clean": True}
    assert {"python", "numpy", "scipy", "torch"} <= set(seal["versions"])
    assert len(seal["eligibility_manifest_sha256"]) == 64


def test_git_probe_unknown_when_git_missing(monkeypatch):
    import subprocess

    def boom(*a, **k):
        raise FileNotFoundError("git")
    monkeypatch.setattr(subprocess, "run", boom)
    monkeypatch.undo()  # drop the autouse probe stub, then break subprocess again
    monkeypatch.setattr(subprocess, "run", boom)
    assert g._git_probe(".", ["x"]) == {"head": "unknown", "clean": "unknown"}


def test_unknown_or_dirty_git_refused_at_entry(setup, monkeypatch):
    cases, folds, _, path = setup
    for state in ({"head": "unknown", "clean": "unknown"}, {"head": "a" * 40, "clean": False}):
        monkeypatch.setattr(g, "_git_probe", lambda root, files, st=state: st)
        with pytest.raises(g.SealMismatch, match="git HEAD unknown or sealed code/protocol files not clean"):
            _call(cases, folds, path)


def test_changed_git_commit_refused_before_any_target_read(setup, monkeypatch):
    _, folds, _, path = setup
    real, reads = _real(_spread_cases())
    monkeypatch.setattr(g, "_git_probe", lambda root, files: {"head": "b" * 40, "clean": True})
    with pytest.raises(g.SealMismatch, match="git commit b{40} differs from sealed a{40}"):
        _call(real, folds, path, allow_synthetic=False, gates=ALL_GATES)
    assert reads == []


def test_seal_built_on_dirty_tree_refused_at_entry(setup):
    cases, folds, _, path = setup
    _tamper(path, lambda s: s.update(git={"head": "a" * 40, "clean": False}))
    with pytest.raises(g.SealMismatch, match="seal: git"):
        _call(cases, folds, path)


def test_expected_seal_sha256_required_and_checked(setup):
    cases, folds, _, path = setup
    with pytest.raises(g.SealMismatch, match="SHA-256"):
        _call(cases, folds, path, expected_seal_sha256="1" * 64)
    with pytest.raises(TypeError):
        g.run_sealed_protocol(cases, folds, seal_path=path, provider_factory=_factory, recipe=TINY,
                              n_bootstrap=N_BOOT, allow_synthetic=True)
    recorded = g.file_sha256(path)
    _tamper(path, lambda s: s.update(created_utc="x"))  # content changed after the hash was recorded
    with pytest.raises(g.SealMismatch, match="SHA-256"):
        _call(cases, folds, path, expected_seal_sha256=recorded)


def test_package_version_mismatch_refused(setup, monkeypatch):
    cases, folds, _, path = setup
    real = g._env_versions()
    for key in ("numpy", "scipy", "torch", "python"):
        monkeypatch.setattr(g, "_env_versions", lambda key=key: {**real, key: "0.0.0-other"})
        with pytest.raises(g.SealMismatch, match="package versions differ"):
            _call(cases, folds, path)


def test_n_bootstrap_below_one_rejected_at_build_and_entry(setup):
    cases, folds, _, path = setup
    for bad in (0, -5):
        with pytest.raises(ValueError, match="n_bootstrap"):
            g.build_seal(cases=cases, fold_by_patient=folds, provider_factory=_factory, recipe=TINY, n_bootstrap=bad)
        with pytest.raises(g.SealMismatch, match="n_bootstrap"):
            _call(cases, folds, path, n_bootstrap=bad)


def test_eligibility_qc_manifest_mismatch_refused(setup):
    cases, folds, _, path = setup
    flipped = [replace(cases[0], qc_passed=not cases[0].qc_passed)] + list(cases[1:])
    with pytest.raises(g.SealMismatch, match="eligibility/QC manifest"):
        _call(flipped, folds, path)
    shifted = [replace(cases[0], forecast_weeks=cases[0].forecast_weeks + 1)] + list(cases[1:])
    with pytest.raises(g.SealMismatch, match="eligibility/QC manifest"):
        _call(shifted, folds, path)


def test_eligibility_manifest_ignores_order_and_voxels():
    cs = _spread_cases()
    a = g.eligibility_manifest_sha256(cs)
    assert a == g.eligibility_manifest_sha256(list(reversed(cs)))
    assert a == g.eligibility_manifest_sha256([replace(c, visit1=c.visit1 * 0) for c in cs])


# ---------------------------------------------------------------- 2026-11-10 guard fixes
def test_seal_records_gates_and_real_run_requires_same_gates(tmp_path):
    cases = _spread_cases()
    real, reads = _real(cases)
    folds = _folds(real)
    seal = g.build_seal(cases=real, fold_by_patient=folds, provider_factory=_factory, recipe=TINY,
                        n_bootstrap=N_BOOT, gates=ALL_GATES)
    assert seal["gates"] == {k: True for k in sorted(ALL_GATES)}
    path = g.write_seal(seal, tmp_path / "seal.json")
    # entry gates differ from sealed gates (an extra gate) -> refused before any target read
    with pytest.raises(g.SealMismatch, match="differ from the sealed gates"):
        _call(real, folds, path, allow_synthetic=False, gates=dict(ALL_GATES, extra_gate=True))
    assert reads == []


def test_seal_without_passed_gates_refused_for_real_cases(setup):
    _, folds, _, path = setup          # seal built with no gates
    real, reads = _real(_spread_cases())
    with pytest.raises(g.SealMismatch, match="sealed gates|seal does not record"):
        _call(real, folds, path, allow_synthetic=False, gates=ALL_GATES)
    assert reads == []


def test_real_case_without_qc_record_refused(setup):
    _, folds, _, path = setup
    real, reads = _real(_spread_cases())
    real[0] = replace(real[0], qc_record=None)
    with pytest.raises(PermissionError, match="no registration QC record"):
        _call(real, folds, path, allow_synthetic=False, gates=ALL_GATES)
    assert reads == []


def test_qc_flag_disagreeing_with_record_refused(setup):
    _, folds, _, path = setup
    real, reads = _real(_spread_cases())
    real[0] = replace(real[0], qc_record=_qc(not real[0].qc_passed))
    with pytest.raises(PermissionError, match="disagrees with its QC record"):
        _call(real, folds, path, allow_synthetic=False, gates=ALL_GATES)
    assert reads == []


def test_qc_record_is_bound_into_manifest():
    real, _ = _real(_spread_cases())
    a = g.eligibility_manifest_sha256(real)
    changed = [replace(real[0], qc_record=dict(real[0].qc_record, reasons=["edited"]))] + real[1:]
    assert a != g.eligibility_manifest_sha256(changed)
    assert g.eligibility_manifest_sha256(_spread_cases()) == g.eligibility_manifest_sha256(list(_spread_cases()))
