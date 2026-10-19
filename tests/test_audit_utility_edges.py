"""Small synthetic edge cases from audit chunk 069; no patient data or output writes."""
import pandas as pd
import pytest

from src.patient_stats import bootstrap_ci, cluster_bootstrap
from src.radiation_model import RadiationSchedule, standard_stupp_schedule


def test_empty_bootstrap_rejects_missing_patient_values():
    with pytest.raises(ValueError, match="at least one patient value"):
        bootstrap_ci([], n_boot=10)


def test_cluster_bootstrap_rejects_missing_patient_groups():
    for frame in (
        pd.DataFrame(columns=["patient_id", "value"]),
        pd.DataFrame({"patient_id": [None], "value": [1.0]}),
    ):
        with pytest.raises(ValueError, match="at least one patient group"):
            cluster_bootstrap(frame, lambda sample: sample["value"].mean(), n_boot=10)


def test_total_radiation_dose_is_invariant_to_fractional_time_origin():
    for start in (0.0, 0.25, 0.5, -0.5):
        assert standard_stupp_schedule(start).total_dose() == 60.0
        assert RadiationSchedule(start, start).total_dose() == 2.0
        assert RadiationSchedule(start, start + 4.9).total_dose() == 10.0


def test_tensor_builder_rejects_nonfinite_inputs_and_regularizes_negative_eigenvalues():
    import numpy as np
    from src.load_ucsf_tensor import build_pd_tensor

    for values in ((float("nan"), 1.0, 0.0), (1.0, float("inf"), 0.0), (1.0, 1.0, float("nan"))):
        with pytest.raises(ValueError, match="must be finite"):
            build_pd_tensor(*values)
    tensor = build_pd_tensor(-2.0, 0.0, 0.6)
    assert np.allclose(tensor, tensor.T)
    assert np.linalg.eigvalsh(tensor).min() >= 1e-6 - 1e-12


def test_report_rankings_follow_supplied_scores_without_reading_patient_artifacts(tmp_path):
    import ast
    from pathlib import Path
    import time
    # Isolate the report function: importing the legacy script creates output/.
    source = Path(__file__).resolve().parents[1] / "src/34_drug_gating_report.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    report = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "generate_report")
    scope = {"Dict": dict, "Path": Path, "time": time}
    exec(compile(ast.Module(body=[report], type_ignores=[]), str(source), "exec"), scope)
    data = {
        "single_ko": [{"gene": "LOW", "collapse_score": 0.01}, {"gene": "HIGH", "collapse_score": 0.9}],
        "single_ko_ti": [{"gene": "SELECTIVE", "therapeutic_index": 3.0, "tumor_collapse": 0.4, "healthy_collapse": 0.0}],
        "dual_ko": [{"gene_a": "A", "gene_b": "B", "bliss_synergy": -0.2}],
        "dual_ko_ti": [{"gene_a": "A", "gene_b": "B", "therapeutic_index": 0.5, "tumor_collapse": 0.7, "healthy_collapse": 0.0, "bliss_synergy": -0.2}],
        "master_switches": [],
    }
    target = tmp_path / "report.md"
    scope["generate_report"](data, str(target))
    rendered = target.read_text(encoding="utf-8")
    assert "HIGH (collapse = 0.9)" in rendered
    assert "SELECTIVE (model score = 3)" in rendered
    assert "A + B (Bliss = -0.2)" in rendered
    assert "TMLHE" not in rendered and "15 pairs" not in rendered
    assert "chosen model cutoffs" in rendered
