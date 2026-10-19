"""Synthetic safety regressions for audit chunk 066; no project pipelines run."""
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent


def load_source(name):
    spec = importlib.util.spec_from_file_location("audit_" + name, ROOT / "src" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_import_graph_builder_does_not_load_latent_arrays(monkeypatch):
    def reject_load(*args, **kwargs):
        raise AssertionError("import attempted graph input I/O")

    monkeypatch.setattr(np, "load", reject_load)
    module = load_source("build_full_graph")
    assert callable(module.main)


def test_import_temporal_comparison_does_not_read_clinical_spreadsheet(monkeypatch):
    def reject_read(*args, **kwargs):
        raise AssertionError("import attempted clinical input I/O")

    monkeypatch.setattr(pd, "read_excel", reject_read)
    module = load_source("compare_naive_vs_treatment")
    assert callable(module.main)


def test_cli_rejects_bulk_execution_and_excludes_audited_unsafe_entries(monkeypatch):
    from src import cli

    assert not set(cli.SCRIPTS.values()).intersection(cli.EXCLUDED_SCRIPTS)
    monkeypatch.setattr(sys, "argv", ["cli", "--run-all"])
    with pytest.raises(SystemExit) as result:
        cli.main()
    assert result.value.code == 2


@pytest.mark.parametrize("name,loader", [
    ("create_4d_3d_embedded", "load_patient_4d_data"),
    ("create_4d_3d_light", "load_patient_4d_data"),
    ("create_4d_3d_viewer", "load_patient_data"),
])
def test_viewers_ignore_cwd_and_retain_explicit_assumed_fallbacks(tmp_path, monkeypatch, name, loader):
    monkeypatch.chdir(tmp_path)
    module = load_source(name)
    assert module.BASE_OUTPUT == ROOT / "output" / "batch_digital_twins"
    patient = tmp_path / "synthetic_patient"
    patient.mkdir()
    np.savez(patient / "clinical_simulation_data.npz", volume_history=np.array([1.0, 2.0]),
             drug_concentration=np.array([0.0, 0.0]), radiation_dose=np.array([0.0, 0.0]))
    (patient / "treatment_schedule.json").write_text(json.dumps({"events": []}), encoding="utf-8")
    monkeypatch.setattr(module, "BASE_OUTPUT", tmp_path)
    data = getattr(module, loader)("synthetic_patient")
    assert data["rho"] == 0.015 and data["D"] == 0.01
    assert data["volume_history"] == [1.0, 2.0]
