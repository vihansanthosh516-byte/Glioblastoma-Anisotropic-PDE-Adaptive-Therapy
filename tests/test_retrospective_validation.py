import numpy as np

from src.retrospective_validation import dice_score, volume_metrics


def test_volume_metrics():
    metrics = volume_metrics([1, 2, 3], [1, 3, 5])
    assert metrics["mae_mm3"] == 1.0
    assert metrics["rmse_mm3"] == np.sqrt(5 / 3)
    assert 0 < metrics["pearson_r"] <= 1


def test_dice_score():
    a = np.array([1, 1, 0, 0])
    b = np.array([1, 0, 1, 0])
    assert dice_score(a, b) == 0.5
