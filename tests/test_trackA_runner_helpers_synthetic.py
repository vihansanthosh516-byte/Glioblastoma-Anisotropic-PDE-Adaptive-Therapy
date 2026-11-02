"""Synthetic-only helper checks for sealed Track A runners; no cohort inputs."""
from importlib import import_module

import numpy as np
import pandas as pd
import pytest


s173 = import_module("src.173_trackA_senmayo_edge_core")
s174 = import_module("src.174_trackA_oxphos_edge_core")
s185 = import_module("src.185_trackA_tam_5ala_core_inv")
s189 = import_module("src.189_trackA_yu2020_lodo")


@pytest.mark.parametrize("module", [s173, s174.s173, s185.s173], ids=["173", "174", "185"])
def test_shared_sign_flip_and_scale_helpers(module):
    p_value, mean = module.sign_flip(np.array([1.0, 1.0]))
    assert p_value == pytest.approx(0.5)
    assert mean == pytest.approx(1.0)

    scaled = module.scale(np.array([3.0, 4.0]))
    assert np.sqrt(np.mean(scaled**2)) == pytest.approx(1.0)


def test_tam_loader_reads_only_synthetic_csvs(tmp_path, monkeypatch):
    sample_ids = ["library P1_core", "library P1_inv", "library P2_core", "library P2_inv"]
    pd.DataFrame(
        {"Annotation": ["Core", "Inv", "Core", "Inv"]}, index=sample_ids
    ).to_csv(tmp_path / "Annotation_col.csv")
    counts = pd.DataFrame(
        {
            sample_ids[0]: [10, 30],
            sample_ids[1]: [20, 20],
            sample_ids[2]: [30, 10],
            sample_ids[3]: [40, 40],
            "width": [100, 100],
        },
        index=["101", "202"],
    )
    counts.to_csv(tmp_path / "Final_count_exon.csv")
    monkeypatch.setattr(s185, "DATA", tmp_path)

    expression, metadata = s185.load()

    assert expression.shape == (2, 4)
    assert metadata.patient.tolist() == ["P1", "P1", "P2", "P2"]
    assert metadata.zone.tolist() == ["core", "inv", "core", "inv"]
    expected = np.log2(np.array([10.0, 30.0]) / 40.0 * 1e6 + 1.0)
    np.testing.assert_allclose(expression["library P1_core"].to_numpy(), expected)


def test_yu_aurc_matches_hand_calculation():
    y = np.array([1, 0, 0, 1])
    p_core = np.array([0.9, 0.8, 0.2, 0.1])
    probabilities = np.column_stack((1.0 - p_core, p_core))
    tie = np.arange(len(y))

    assert s189.class_balanced_aurc(y, probabilities, tie) == pytest.approx(0.34375)
    record = s189.metric_record(y, probabilities, tie)
    assert record["balanced_accuracy"] == pytest.approx(0.5)
    assert record["brier"] == pytest.approx(0.375)


def test_yu_probability_and_tuning_guards_reject_bad_inputs():
    with pytest.raises(s189.StopRun, match="sum to one"):
        s189.checked_probabilities(np.array([[0.2, 0.2]]), 1)

    scores = np.array([[1.0 + 5e-13, 1.0], [1.0 + 4e-13, 3.0]])
    assert s189.first_setting_within_tolerance(scores) == (0, 0)


def test_yu_tie_ranks_are_deterministic_and_normalize_zero_rows():
    ids = np.array(["cell-b", "cell-a", "cell-c"])
    first = s189.tie_ranks(ids)
    second = s189.tie_ranks(ids)
    np.testing.assert_array_equal(first, second)
    np.testing.assert_array_equal(np.sort(first), np.arange(len(ids)))

    normalized = s189.normalize_rows(np.array([[3.0, 4.0], [0.0, 0.0]]))
    np.testing.assert_allclose(normalized[0], [0.6, 0.8])
    np.testing.assert_array_equal(normalized[1], [0.0, 0.0])
