"""Track A external-test lock machinery (handoff section 5). No real data here.

Stages: preflight -> predict -> sealed evaluate. Each stage checks the signed
protocol, config and input hashes. Predictions are saved before any outcome
file is read. Outputs are written atomically and never overwritten.
"""

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats


# ------------------------------------------------------------- scoring units
def regional_unit_scores(margins, labels, *, tie_tol=1e-9):
    """Per patient: 1 correct paired order, 0.5 tie, 0 wrong.

    margins: {patient: signed margin, positive = infiltrating ranked above core}.
    labels: {patient: +1 or -1, the pathology order}.
    """
    if set(margins) != set(labels):
        raise ValueError("Every enrolled patient needs a prediction, including fallback cases")
    out = {}
    for p, m in margins.items():
        if not np.isfinite(m) or labels[p] not in (-1, 1):
            raise ValueError(f"Invalid margin or label for {p}")
        s = m * labels[p]
        out[p] = 0.5 if abs(m) <= tie_tol else (1.0 if s > 0 else 0.0)
    return out


def tmz_unit_errors(predicted, observed):
    """Per origin absolute error of the origin-mean viability."""
    if set(predicted) != set(observed):
        raise ValueError("Every enrolled origin needs a prediction, including fallback cases")
    return {o: abs(float(predicted[o]) - float(observed[o])) for o in predicted}


def paired_differences(endpoint, candidate, comparator):
    """Larger favours the candidate for both endpoints."""
    if set(candidate) != set(comparator):
        raise ValueError("Arms must be scored on the same units")
    keys = sorted(candidate)
    if endpoint == "regional":
        return np.array([candidate[k] - comparator[k] for k in keys])
    if endpoint == "tmz":
        return np.array([comparator[k] - candidate[k] for k in keys])
    raise ValueError("endpoint must be 'regional' or 'tmz'")


# --------------------------------------------------------------- statistics
def lower_bound(diffs, alpha):
    """One-sided t lower bound for the mean paired improvement."""
    d = np.asarray(diffs, float)
    if d.size < 2:
        return -np.inf
    sd = d.std(ddof=1)
    if sd == 0:
        return float(d.mean())
    return float(d.mean() - stats.t.ppf(1 - alpha, d.size - 1) * sd / np.sqrt(d.size))


def sign_test(diffs):
    """Exact one-sided sign test on discordant units (wins vs losses only)."""
    d = np.asarray(diffs, float)
    wins, losses = int((d > 0).sum()), int((d < 0).sum())
    n = wins + losses
    p = 1.0 if n == 0 else float(stats.binomtest(wins, n, 0.5, alternative="greater").pvalue)
    return {"wins": wins, "losses": losses, "ties": int(d.size - n), "p": p}


def holm(pvalues):
    names = sorted(pvalues, key=pvalues.get)
    adjusted, running = {}, 0.0
    for i, name in enumerate(names):
        running = max(running, min(1.0, (len(names) - i) * pvalues[name]))
        adjusted[name] = running
    return adjusted


def decide(diffs_by_comparator, *, min_gain, alpha):
    """Pass only if every primary comparator's Bonferroni lower bound exceeds min_gain."""
    if not diffs_by_comparator:
        raise ValueError("A primary comparator family is required")
    level = alpha / len(diffs_by_comparator)
    detail = {}
    for name, d in diffs_by_comparator.items():
        d = np.asarray(d, float)
        detail[name] = {"n": int(d.size), "mean": float(d.mean()), "lower_bound": lower_bound(d, level),
                        "sign_test": sign_test(d)}
    adj = holm({k: v["sign_test"]["p"] for k, v in detail.items()})
    for k in detail:
        detail[k]["sign_test"]["holm_p"] = adj[k]
    passed = all(v["lower_bound"] > min_gain for v in detail.values())
    return {"passed": bool(passed), "min_gain": min_gain, "alpha": alpha, "per_comparator_alpha": level,
            "comparators": detail}


def power(n_units, *, mean_gain, sd, n_comparators, min_gain, alpha, correlation=0.5,
          tie_probability=0.0, attrition=0.0, reps=2000, seed=0):
    """Monte Carlo probability that `decide` passes under a correlated normal paired model."""
    rng = np.random.default_rng(seed)
    cov = sd**2 * (correlation * np.ones((n_comparators, n_comparators))
                   + (1 - correlation) * np.eye(n_comparators))
    hits = 0
    for _ in range(reps):
        n = rng.binomial(n_units, 1 - attrition)
        if n < 2:
            continue
        d = rng.multivariate_normal(np.full(n_comparators, mean_gain), cov, size=n)
        d[rng.random(n) < tie_probability] = 0.0
        hits += decide({f"c{j}": d[:, j] for j in range(n_comparators)}, min_gain=min_gain, alpha=alpha)["passed"]
    return hits / reps


def units_for_power(target, *, start=4, stop=2000, **kw):
    """Smallest unit count reaching `target` power, or None."""
    for n in range(start, stop + 1):
        if power(n, **kw) >= target:
            return n
    return None


# -------------------------------------------------------------- lock + files
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def atomic_write_json(path, payload):
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"{path} exists; reruns are not allowed")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    with open(tmp, "x", encoding="utf-8") as f:
        json.dump(payload, f, indent=1, allow_nan=False)
        f.flush()
        os.fsync(f.fileno())
    os.rename(tmp, path)  # fails if path appeared meanwhile on Windows


REQUIRED_PROTOCOL = ("protocol_id", "endpoint", "signed", "signer", "commit", "config_sha256",
                     "inputs", "outcomes", "primary_comparators", "gates")


def verify_lock(protocol_path, config_path):
    """Return the protocol if signed and every hash matches; raise otherwise."""
    protocol = json.loads(Path(protocol_path).read_text(encoding="utf-8"))
    missing = [k for k in REQUIRED_PROTOCOL if k not in protocol]
    if missing:
        raise ValueError(f"Protocol missing fields: {missing}")
    if protocol["signed"] is not True or not protocol["signer"]:
        raise PermissionError("Protocol is not signed off")
    if sha256(config_path) != protocol["config_sha256"]:
        raise ValueError("Config hash does not match the signed protocol")
    for path, digest in protocol["inputs"].items():
        if sha256(path) != digest:
            raise ValueError(f"Input hash changed: {path}")
    for key in ("min_gain", "alpha"):
        if key not in protocol["gates"] or protocol["gates"][key] is None:
            raise ValueError(f"Gate {key} must be fixed before unsealing")
    return protocol


def preflight(manifest_rows, *, origin_key="canonical_origin", role_key="split_role"):
    """Outcome-free checks: one role per origin, no unknown roles, eval roles disjoint."""
    roles = {}
    for row in manifest_rows:
        o, r = row[origin_key], row[role_key]
        if r not in ("fit", "calibration", "primary_external", "shift_external", "excluded"):
            raise ValueError(f"Unknown split role {r!r}")
        if roles.setdefault(o, r) != r:
            raise ValueError(f"Origin {o} appears in two roles")
    return {r: sorted(o for o, x in roles.items() if x == r) for r in set(roles.values())}


def save_predictions(run_dir, protocol, predictions):
    """predictions: {arm: {unit: value}}. Stored before outcomes are read."""
    for arm, values in predictions.items():
        if any(k in ("viability", "label", "outcome") for k in values):
            raise ValueError("Outcome-like keys in predictions")
        if not all(np.isfinite(float(v)) for v in values.values()):
            raise ValueError(f"Non-finite prediction in arm {arm}")
    units = {frozenset(v) for v in predictions.values()}
    if len(units) != 1:
        raise ValueError("All arms must predict the same units")
    path = Path(run_dir) / "predictions.json"
    atomic_write_json(path, {"protocol_id": protocol["protocol_id"], "predictions": predictions})
    return sha256(path)


def sealed_evaluate(run_dir, protocol_path, config_path, candidate):
    """Verify lock, require saved predictions, log outcome access, score once."""
    run_dir = Path(run_dir)
    protocol = verify_lock(protocol_path, config_path)
    pred_path = run_dir / "predictions.json"
    if not pred_path.exists():
        raise RuntimeError("Predictions must be saved before outcomes are opened")
    if (run_dir / "result.json").exists() or (run_dir / "outcome_access.json").exists():
        raise FileExistsError("This sealed evaluation already ran; reruns are not allowed")
    saved = json.loads(pred_path.read_text(encoding="utf-8"))
    if saved["protocol_id"] != protocol["protocol_id"]:
        raise ValueError("Predictions belong to a different protocol")
    atomic_write_json(run_dir / "outcome_access.json", {
        "protocol_id": protocol["protocol_id"], "predictions_sha256": sha256(pred_path),
        "opened_at": datetime.now(timezone.utc).isoformat(), "outcomes": protocol["outcomes"]})
    outcomes = json.loads(Path(protocol["outcomes"]).read_text(encoding="utf-8"))
    preds, endpoint = saved["predictions"], protocol["endpoint"]
    score = (lambda arm: regional_unit_scores(preds[arm], {k: int(v) for k, v in outcomes.items()})) \
        if endpoint == "regional" else (lambda arm: tmz_unit_errors(preds[arm], outcomes))
    cand = score(candidate)
    diffs = {c: paired_differences(endpoint, cand, score(c)) for c in protocol["primary_comparators"]}
    decision = decide(diffs, min_gain=protocol["gates"]["min_gain"], alpha=protocol["gates"]["alpha"])
    result = {"protocol_id": protocol["protocol_id"], "candidate": candidate, "units": len(cand),
              "decision": decision,
              "paired_differences": {c: d.tolist() for c, d in diffs.items()}}
    atomic_write_json(run_dir / "result.json", result)
    return result
