import numpy as np
import pytest

from trackc_world_summary import paired_world_summary


def test_equal_world_summary_and_seeded_interval():
    candidate = np.array([[1.0, 3.0], [10.0, 14.0], [5.0, 7.0]])
    comparator = np.array([[0.0, 0.0], [8.0, 8.0], [5.0, 5.0]])
    args = (candidate, comparator, ["a", "b", "c"])
    first = paired_world_summary(*args, resamples=2000, seed=11)
    second = paired_world_summary(*args, resamples=2000, seed=11)
    np.testing.assert_allclose(first["world_mean_differences"], [2, 4, 1])
    assert first["equal_world_mean_difference"] == pytest.approx(7 / 3)
    np.testing.assert_array_equal(first["bootstrap_percentile_95"], second["bootstrap_percentile_95"])
    assert first["bootstrap_unit"] == "world"


@pytest.mark.parametrize("candidate,comparator,ids", [
    ([[1]], [[0]], ["a"]),
    ([[1], [2]], [[0]], ["a", "b"]),
    ([[1], [np.nan]], [[0], [0]], ["a", "b"]),
    ([[1], [2]], [[0], [0]], ["a", "a"]),
])
def test_invalid_pairs_are_rejected(candidate, comparator, ids):
    with pytest.raises(ValueError):
        paired_world_summary(candidate, comparator, ids, resamples=10, seed=1)
