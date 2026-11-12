import math
import unittest

from track_b_acrin6677_forecast import fit_m, predict_b2, predict_m, score
from track_b_acrin6677_run import build_windows, exact_wilcoxon, patient_errors


class ForecastTests(unittest.TestCase):
    def test_model_recovers_both_branches(self):
        for g, k in [(0.02, 0), (0, 0.03)]:
            def curve(t):
                return 2 + g * t + k * 56 * math.expm1(-t / 56)

            fitted = fit_m(5, curve(0), 25, curve(20))
            self.assertAlmostEqual(fitted[0], g)
            self.assertAlmostEqual(fitted[1], k)
            self.assertAlmostEqual(predict_m(5, curve(0), 25, curve(20), 55), curve(50))

    def test_rising_matches_b2(self):
        args = (0, 1, 5, 2, 20)
        self.assertEqual(predict_m(*args), predict_b2(*args))

    def test_finite(self):
        for y2 in (0, 1, 2):
            self.assertTrue(math.isfinite(predict_m(0, 1, 5, y2, 100)))

    def test_three_patient_scores(self):
        for pred, obs, expected in [(4, 1, 3), (1, 3, 2), (2, 2, 0)]:
            self.assertEqual(score(pred, obs), expected)
        flat = ((0, 4), (1, 4), (2, 4))
        self.assertEqual(list(patient_errors(flat).values()), [0, 0, 0, 0])

    def test_windows(self):
        def row(patient, day, area, no="", quality=""):
            return {"cn": patient, "S0e2d": str(day), "s0e23": str(area),
                    "s0e13": no, "s0e5": quality}
        rows = [row("a", 5, 2), row("a", 3, 4), row("a", 3, 6),
                row("a", 1, "", "1"), row("a", 2, 20, quality="1"),
                row("a", 8, 50), row("b", 0, 1), row("b", 4, 1)]
        self.assertEqual(build_windows(rows), {"a": ((1, 0), (3, 5), (5, 2))})

    def test_exact_signed_ranks_with_ties_and_zeros(self):
        self.assertEqual(exact_wilcoxon([0, 0]), (1, False))
        self.assertEqual(exact_wilcoxon([0, 1, 1, 2]), (0.25, True))
        self.assertEqual(exact_wilcoxon([0, -1, 1, 2]), (0.75, True))


if __name__ == "__main__":
    unittest.main()
