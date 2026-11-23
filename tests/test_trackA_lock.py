import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from src import trackA_lock as L


def _setup(tmp, endpoint="tmz", signed=True, gates=None):
    tmp = Path(tmp)
    cfg = tmp / "config.json"
    cfg.write_text('{"penalty": 0.2}')
    inp = tmp / "expr.csv"
    inp.write_text("origin,g1\nA,1\nB,2\nC,3\n")
    out = tmp / "outcomes.json"
    if endpoint == "tmz":
        out.write_text(json.dumps({"A": 0.5, "B": 0.6, "C": 0.7}))
    else:
        out.write_text(json.dumps({"A": 1, "B": -1, "C": 1}))
    proto = tmp / "protocol.json"
    proto.write_text(json.dumps({
        "protocol_id": "P1", "endpoint": endpoint, "signed": signed, "signer": "reviewer",
        "commit": "abc", "config_sha256": L.sha256(cfg), "inputs": {str(inp): L.sha256(inp)},
        "outcomes": str(out), "primary_comparators": ["ridge"],
        "gates": gates or {"min_gain": 0.0, "alpha": 0.05}}))
    return proto, cfg, inp, tmp / "run"


class Scoring(unittest.TestCase):
    def test_regional_ties_and_missing(self):
        s = L.regional_unit_scores({"a": 1.0, "b": 0.0, "c": -2.0}, {"a": 1, "b": 1, "c": 1})
        self.assertEqual(s, {"a": 1.0, "b": 0.5, "c": 0.0})
        with self.assertRaises(ValueError):
            L.regional_unit_scores({"a": 1.0}, {"a": 1, "b": 1})

    def test_tmz_direction(self):
        cand = L.tmz_unit_errors({"a": 0.5}, {"a": 0.6})
        comp = L.tmz_unit_errors({"a": 0.9}, {"a": 0.6})
        self.assertAlmostEqual(L.paired_differences("tmz", cand, comp)[0], 0.2)

    def test_holm_monotone(self):
        adj = L.holm({"x": 0.01, "y": 0.04, "z": 0.03})
        self.assertAlmostEqual(adj["x"], 0.03)
        self.assertAlmostEqual(adj["z"], 0.06)
        self.assertAlmostEqual(adj["y"], 0.06)

    def test_decide_needs_every_comparator(self):
        good, bad = np.full(30, 0.2) + np.linspace(-0.01, 0.01, 30), np.linspace(-0.1, 0.1, 30)
        self.assertTrue(L.decide({"a": good}, min_gain=0.1, alpha=0.05)["passed"])
        self.assertFalse(L.decide({"a": good, "b": bad}, min_gain=0.1, alpha=0.05)["passed"])

    def test_power_grows_with_n(self):
        kw = dict(mean_gain=0.05, sd=0.1, n_comparators=2, min_gain=0.0, alpha=0.05, reps=300, seed=1)
        self.assertLess(L.power(10, **kw), L.power(80, **kw))
        self.assertLess(L.power(80, **{**kw, "mean_gain": 0.0}), 0.1)


class Lock(unittest.TestCase):
    def test_full_run_and_no_rerun(self):
        with tempfile.TemporaryDirectory() as t:
            proto, cfg, _, run = _setup(t)
            L.save_predictions(run, {"protocol_id": "P1"},
                               {"cand": {"A": 0.5, "B": 0.6, "C": 0.7}, "ridge": {"A": 0.9, "B": 0.2, "C": 0.3}})
            res = L.sealed_evaluate(run, proto, cfg, "cand")
            self.assertEqual(res["units"], 3)
            self.assertTrue((run / "outcome_access.json").exists())
            with self.assertRaises(FileExistsError):
                L.sealed_evaluate(run, proto, cfg, "cand")

    def test_unsigned_changed_hash_missing_gate(self):
        with tempfile.TemporaryDirectory() as t:
            proto, cfg, _, _ = _setup(t, signed=False)
            with self.assertRaises(PermissionError):
                L.verify_lock(proto, cfg)
        with tempfile.TemporaryDirectory() as t:
            proto, cfg, inp, _ = _setup(t)
            inp.write_text("changed")
            with self.assertRaises(ValueError):
                L.verify_lock(proto, cfg)
        with tempfile.TemporaryDirectory() as t:
            proto, cfg, _, _ = _setup(t, gates={"min_gain": None, "alpha": 0.05})
            with self.assertRaises(ValueError):
                L.verify_lock(proto, cfg)

    def test_outcomes_need_saved_predictions(self):
        with tempfile.TemporaryDirectory() as t:
            proto, cfg, _, run = _setup(t)
            with self.assertRaises(RuntimeError):
                L.sealed_evaluate(run, proto, cfg, "cand")
            self.assertFalse((run / "outcome_access.json").exists())

    def test_prediction_guards(self):
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaises(ValueError):
                L.save_predictions(t, {"protocol_id": "P"}, {"a": {"x": 1.0}, "b": {"y": 1.0}})
            with self.assertRaises(ValueError):
                L.save_predictions(t, {"protocol_id": "P"}, {"a": {"viability": 1.0}})
            with self.assertRaises(ValueError):
                L.save_predictions(t, {"protocol_id": "P"}, {"a": {"x": float("nan")}})

    def test_interrupted_write_leaves_no_result(self):
        with tempfile.TemporaryDirectory() as t:
            target = Path(t) / "r.json"
            with self.assertRaises(ValueError):
                L.atomic_write_json(target, {"x": float("nan")})
            self.assertFalse(target.exists())

    def test_preflight_duplicate_origin(self):
        rows = [{"canonical_origin": "A", "split_role": "fit"}, {"canonical_origin": "A", "split_role": "primary_external"}]
        with self.assertRaises(ValueError):
            L.preflight(rows)
        ok = L.preflight([{"canonical_origin": "A", "split_role": "fit"}, {"canonical_origin": "B", "split_role": "primary_external"}])
        self.assertEqual(ok["primary_external"], ["B"])


if __name__ == "__main__":
    unittest.main()
