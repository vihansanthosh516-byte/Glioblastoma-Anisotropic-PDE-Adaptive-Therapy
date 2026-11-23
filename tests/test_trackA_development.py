"""Regression checks for the Track A development paths, using generated inputs."""

import unittest

import numpy as np
from sklearn.linear_model import Ridge

from src.trackA_cohort_preflight import audit
from src.trackA_development_benchmark import _regional_fixture, run
from src.trackA_regional_model import calibrate_margin_bound, fit_regional, spatial_fallback
from src.trackA_tmz_model import (
    Cultures, _fit_shared, benchmark_tmz, matched_signature_weights, predict_tmz_arms, score_origins,
)


class TestRegionalDevelopment(unittest.TestCase):
    def test_fallback_reasons_keep_every_pair(self):
        fallback = np.array([-5.0] * 8)
        margins = [[2, 3], [-2, -3], [0.5, 2], [2, -3], [2, 3], [np.nan, 3], [2, 3], [0, 3]]
        decision = spatial_fallback(
            fallback, margins, bound=[1, 1, 0.5, 0, np.inf, 0, 0, 0],
            audit_qualified=[True, True, True, True, True, True, False, True],
        )
        np.testing.assert_array_equal(decision["scores"]["bounded_spatial"], [2, -2, -5, -5, -5, -5, -5, -5])
        np.testing.assert_array_equal(decision["used_spatial"]["bounded_spatial"], [True, True] + [False] * 6)
        self.assertEqual(decision["scores"]["agreement_only"][6], 2.0)

    def test_calibration_uses_joint_maximum_and_finite_sample_rank(self):
        margins = [[-0.2, 0.1], [0.15, 0.35], [-0.3, -0.6], [0.8, -0.4]]
        self.assertEqual(calibrate_margin_bound(margins, [0, 0, 0, 0], alpha=0.4), 0.6)
        self.assertTrue(np.isinf(calibrate_margin_bound(margins, [0, 0, 0, 0], alpha=0.1)))

    def test_pair_order_reverses_every_margin_without_changing_coverage(self):
        rng = np.random.default_rng(11)
        features, order, _ = _regional_fixture(rng, 40)
        model = fit_regional(features, order, strength=0.1)
        bound, audit_qualified = 0.25, np.ones(40, dtype=bool)
        scores, used = model.predict(features, bound=bound, audit_qualified=audit_qualified)
        reverse, reverse_used = model.predict(
            {name: -values for name, values in features.items()}, bound=bound, audit_qualified=audit_qualified,
        )
        for arm in scores:
            np.testing.assert_allclose(scores[arm], -reverse[arm])
            np.testing.assert_array_equal(used[arm], reverse_used[arm])


class TestTmzDevelopment(unittest.TestCase):
    def test_joint_loss_matches_weighted_normal_equation(self):
        x = np.array([[0], [2], [4], [0], [6], [7]], dtype=float)
        y = np.array([0, 1, 2, 3, 9, 8], dtype=float)
        groups = np.array(["a"] * 3 + ["b"] * 2 + ["c"])
        penalty, within_weight = 0.2, 0.8
        mean_x = np.array([2.0, 3.0, 7.0])
        mean_y = np.array([1.0, 6.0, 8.0])
        between_x, between_y = mean_x - mean_x.mean(), mean_y - mean_y.mean()
        variance = np.mean(between_x**2) + within_weight * (8 / 3 + 9) / 2
        covariance = np.mean(between_x * between_y) + within_weight * (4 / 3 + 9) / 2
        expected = covariance / (variance + penalty)
        intercept, coefficient = _fit_shared(x, y, groups, groups, penalty, within_weight)
        np.testing.assert_allclose(coefficient, [expected])
        self.assertAlmostEqual(intercept, mean_y.mean() - mean_x.mean() * expected)

    def test_origin_balanced_ridge_equals_a_special_case_of_joint_loss(self):
        x = np.array([[0], [2], [4], [0], [6], [7]], dtype=float)
        y = np.array([0, 1, 2, 3, 9, 8], dtype=float)
        groups = np.array(["a"] * 3 + ["b"] * 2 + ["c"])
        penalty = 0.2
        weights = np.array([1 / 9] * 3 + [1 / 6] * 2 + [1 / 3])
        baseline = Ridge(alpha=penalty).fit(x, y, sample_weight=weights)
        intercept, coefficient = _fit_shared(x, y, groups, groups, penalty, within_weight=2 / 3)
        np.testing.assert_allclose(coefficient, baseline.coef_)
        self.assertAlmostEqual(intercept, baseline.intercept_)

    def test_equal_origin_weight_and_singletons(self):
        data = Cultures([[0], [0], [0], [0]], [0, 2, 4, 10], ["a", "a", "a", "b"], ["1", "2", "3", "4"])
        scored = score_origins(data, {"arm": [2, 2, 2, 2]})["arm"]
        self.assertEqual(scored["origin_mae"], 4.0)
        self.assertAlmostEqual(scored["mean_within_origin_rmse"], np.sqrt(8 / 3))
        self.assertIsNone(scored["per_origin"]["b"]["within_origin_rmse"])

    def test_reject_duplicate_cultures_and_cross_split_origins(self):
        with self.assertRaisesRegex(ValueError, "duplicate cultures"):
            Cultures([[0], [0]], [0, 1], ["a", "a"], ["1", "1"])
        training = Cultures([[0], [1]], [0, 1], ["a", "a"], ["1", "2"])
        evaluation = Cultures([[2]], [2], ["a"], ["3"])
        with self.assertRaisesRegex(ValueError, "disjoint"):
            benchmark_tmz(training, evaluation, penalty=0.2, within_weight=1.0,
                          strata=["same", "same"], shuffle_seed=0)

    def test_signatures_preserve_bin_and_weight_and_exclude_source_genes(self):
        fixed = [1, 0, 0, -2, 0, 0]
        rows = matched_signature_weights(fixed, ["a"] * 3 + ["b"] * 3, count=10, seed=4)
        np.testing.assert_array_equal(rows[:, [0, 3]], np.zeros((10, 2)))
        np.testing.assert_array_equal(rows[:, :3].sum(axis=1), np.ones(10))
        np.testing.assert_array_equal(rows[:, 3:].sum(axis=1), np.full(10, -2))

    def test_blinded_prediction_requires_no_evaluation_viability(self):
        training = Cultures(
            [[0], [2], [4], [0], [6], [7]], [0, 1, 2, 3, 9, 8],
            ["a", "a", "a", "b", "b", "c"], ["1", "2", "3", "4", "5", "6"],
        )
        result = predict_tmz_arms(
            training, [[1], [5]], ["d", "d"], ["7", "8"],
            penalty=0.2, within_weight=1, strata=["same"] * 6, shuffle_seed=4,
        )
        self.assertEqual(result["within_origin"].shape, (2,))
        self.assertTrue(np.isfinite(result["within_origin"]).all())


class TestMetadataPreflight(unittest.TestCase):
    def test_regional_counts_people_and_excludes_expression_derived_labels(self):
        rows = []
        for role, assay in (("development", "a"), ("calibration", "a"), ("evaluation", "a"), ("shift", "b")):
            for region in ("core", "infiltrating"):
                rows.append(dict(
                    study="dev" if role == "calibration" else role.replace("development", "dev"),
                    accession=f"synthetic-{role}", origin=role, item=f"{role}-{region}", role=role,
                    condition=assay, reuse_terms="verified", origin_lineage="verified",
                    prior_exposure="untouched", evidence="synthetic", section="slide",
                    region=region, label_basis="independent_pathology", modality="cell_spatial",
                    registration="verified", spatial_inputs="verified", independent_audit="verified",
                ))
        rows.extend([
            {**rows[4], "item": "evaluation-core-replicate"},
            {**rows[4], "origin": "invalid", "item": "invalid-core", "label_basis": "other"},
        ])
        result = audit(rows, lane="cell-spatial", minimum_evaluation=1, minimum_shift=1,
                       minimum_repeated_development=None, comparison_count=2, alpha=0.05)
        self.assertEqual(result["roles"]["evaluation"]["eligible_origins"], ["evaluation"])
        self.assertEqual(result["excluded_items"][0]["reasons"], ["label_basis"])
        self.assertTrue(result["checks"]["independent_calibration_pairs"])
        self.assertFalse(result["checks"]["sign_test_attainability"])

    def test_tmz_inventory_requires_real_sister_cultures_and_untouched_origins(self):
        base = dict(
            study="dev", accession="synthetic-a", origin="a", item="culture1", role="development",
            condition="plate_A", reuse_terms="verified", origin_lineage="verified",
            prior_exposure="seen", evidence="synthetic", untreated_rna="verified",
            measured_viability="verified", assay_definition="verified", independent_culture="verified",
            mgmt="unknown",
        )
        rows = [
            base,
            {**base, "item": "well2", "independent_culture": "missing"},
            {**base, "study": "eval", "accession": "synthetic-b", "role": "evaluation",
             "origin": "b", "item": "culture2", "prior_exposure": "unknown"},
        ]
        settings = dict(lane="tmz", minimum_evaluation=1, minimum_shift=1, minimum_repeated_development=1)
        result = audit(rows, **settings)
        self.assertFalse(result["checks"]["independent_sister_cultures"])
        self.assertFalse(result["checks"]["evaluation_origin_floor"])
        rows.append({**base, "item": "culture4"})
        with_sisters = audit(rows, **settings)
        self.assertTrue(with_sisters["checks"]["independent_sister_cultures"])
        self.assertEqual(with_sisters["roles"]["development"]["mgmt_complete_origins"], [])
        with self.assertRaisesRegex(ValueError, "crosses a role"):
            audit([*rows, {**rows[2], "origin": "a", "item": "culture5"}], **settings)

    def test_one_study_cannot_obtain_freshness_via_a_second_accession(self):
        base = dict(study="same", accession="first", origin="a", item="a", role="development",
                    condition="r1", reuse_terms="verified", origin_lineage="verified",
                    prior_exposure="untouched", evidence="synthetic", region="core",
                    label_basis="independent_pathology", modality="bulk_rna")
        with self.assertRaisesRegex(ValueError, "crosses study roles"):
            audit([base, {**base, "accession": "second", "origin": "b", "item": "b", "role": "evaluation"}],
                  lane="whole-tissue", minimum_evaluation=1, minimum_shift=1,
                  minimum_repeated_development=None, comparison_count=1, alpha=0.05)


class TestDevelopmentRun(unittest.TestCase):
    def test_all_scenarios_exercise_the_intended_paths(self):
        report = run()
        self.assertTrue(all(report["checks"].values()), report["checks"])
        primary = report["tmz"]["signal"]["within_origin"]["per_origin"]
        shift = report["tmz"]["assay_shift"]["within_origin"]["per_origin"]
        self.assertEqual({k: v["observed"] for k, v in primary.items()},
                         {k: v["observed"] for k, v in shift.items()})


if __name__ == "__main__":
    unittest.main()
