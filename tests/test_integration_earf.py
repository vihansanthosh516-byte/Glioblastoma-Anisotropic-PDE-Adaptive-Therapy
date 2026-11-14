import numpy as np
from trackc_research_grid import research_model_grid
import integration_earf as e


def test_set_is_subset_and_drops_strong_models():
    grid = research_model_grid()
    s = e.evidence_anchored_model_set(grid)
    assert 0 < len(s) < len(grid) and set(s) <= set(grid)
    assert all(abs(m.a) < 0.5 or abs(m.b) < 0.5 for m in s)


def test_weights_normalized_and_deterministic():
    grid = research_model_grid()
    w1, w2 = e.evidence_weights(grid), e.evidence_weights(grid)
    assert np.isclose(w1.sum(), 1) and np.array_equal(w1, w2)


def test_sampler_in_bounds():
    rng = np.random.default_rng(0)
    for _ in range(200):
        m = e.sample_evidence_world(rng, [0.6, 0.7, 0.8, 0.9])
        assert -1 <= m.a <= 1 and -1 <= m.b <= 1 and -0.95 <= m.d <= 0.95


def test_credible_box_size_and_range():
    s = e.credible_box_model_set([0.5, 0.6, 0.7, 0.8, 0.9])
    assert len(s) == 135 and len(set(s)) == 135
    assert all(abs(m.a) < 0.35 and abs(m.b) < 0.4 for m in s)
