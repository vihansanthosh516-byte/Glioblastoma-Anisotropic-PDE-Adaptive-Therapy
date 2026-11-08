import json

import numpy as np
import pytest
from scipy.stats import t

from trackc_t_bonferroni_coverage_audit import (
    DISTRIBUTION_PAIRS,
    WORLD_COUNTS,
    bonferroni_t_intervals,
    draw_standardized,
    wilson_interval,
)


def test_distribution_catalog_and_standardization():
    assert len(DISTRIBUTION_PAIRS) == 6
    assert WORLD_COUNTS == (5, 10, 14, 20)
    rng = np.random.Generator(np.random.PCG64DXSM(4321))
    for name in ("normal", "centered_gamma2", "uniform"):
        values = draw_standardized(name, 100_000, rng)
        assert values.shape == (100_000,)
        assert np.isfinite(values).all()
        assert abs(values.mean()) < 0.025
        assert abs(values.var() - 1.0) < 0.06
    t3_values = draw_standardized("student_t3", 100_000, rng)
    assert abs(np.mean(np.abs(t3_values) <= 1.0 / np.sqrt(3.0)) - 0.609) < 0.02


def test_two_family_intervals_use_bonferroni_student_t_quantile():
    samples = {
        "single-response": np.array([-2.0, -1.0, 1.0, 2.0]),
        "mixture-response": np.array([-3.0, 0.0, 1.0, 2.0]),
    }
    result = bonferroni_t_intervals(samples, confidence_level=0.95)
    expected_critical = t.ppf(1 - 0.05 / (2 * 2), df=3)
    assert result["critical_value"] == pytest.approx(expected_critical)
    for name, values in samples.items():
        row = result["intervals"][name]
        expected_se = values.std(ddof=1) / np.sqrt(values.size)
        assert row["mean"] == pytest.approx(values.mean())
        assert row["standard_error"] == pytest.approx(expected_se)
        assert row["lower"] == pytest.approx(values.mean() - expected_critical * expected_se)
        assert row["upper"] == pytest.approx(values.mean() + expected_critical * expected_se)


@pytest.mark.parametrize("samples", [
    {"single-response": [1.0, 2.0]},
    {"single-response": [1.0, 2.0], "mixture-response": [1.0, 1.0]},
    {"single-response": [1.0, np.nan], "mixture-response": [1.0, 2.0]},
    {"single-response": [1.0, 2.0], "other": [1.0, 2.0]},
])
def test_inference_rejects_missing_bad_or_degenerate_family(samples):
    with pytest.raises(ValueError):
        bonferroni_t_intervals(samples)


def test_wilson_interval_handles_edges_and_midpoint():
    lo, hi = wilson_interval(0, 100)
    assert lo == 0.0
    assert hi > 0.0
    lo, hi = wilson_interval(50, 100)
    assert 0.40 < lo < 0.42
    assert 0.58 < hi < 0.60


def test_main_writes_complete_audit_and_refuses_overwrite(tmp_path, monkeypatch):
    import trackc_t_bonferroni_coverage_audit as audit

    output_path = tmp_path / "coverage.json"
    monkeypatch.setattr(audit, "OUTPUT_PATH", output_path)
    audit.main()
    result = json.loads(output_path.read_text(encoding="utf-8"))
    assert result["record_type"] == "outcome_blind_trackc_t_bonferroni_coverage_screen"
    assert result["cells_total"] == 24
    assert result["inference_failures_total"] == 0
    assert result["patient_data_used"] is False
    assert result["policy_outcomes_used"] is False
    assert result["provenance"]["code_sha256"]
    assert result["provenance"]["protocol_sha256"]
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        audit.main()
