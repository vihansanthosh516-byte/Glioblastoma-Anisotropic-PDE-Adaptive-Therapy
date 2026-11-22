import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from trackA_order_simulation import draw, contrasts, concordance, LogisticRegression


def blocks(x):
    return np.log1p(10000 * x / x.sum(axis=-1, keepdims=True))


def corrupt(rng, x):
    i, j = np.arange(len(x))[:, None], np.arange(2)
    k = rng.integers(4, size=(len(x), 2))
    original = x[i, j, k]
    p = original * np.exp(rng.normal(scale=3, size=original.shape))
    x[i, j, k] = rng.multinomial(original.sum(axis=-1), p / p.sum(axis=-1, keepdims=True))


def run(train, test):
    a, b = blocks(train), blocks(test)
    train_features, test_features = contrasts(train), contrasts(test)
    delta = train_features["block_logistic"]
    model = LogisticRegression(C=1, fit_intercept=False).fit(
        np.concatenate([delta, -delta]), np.repeat([1, 0], len(delta)))

    def centered(x):
        scores = x @ model.coef_[0]
        return scores - np.median(scores, axis=(1, 2), keepdims=True)

    scale = np.sqrt(np.mean(centered(a)**2)) or 1
    s = centered(b) / scale
    pools = {
        "bounded": logsumexp(np.clip(s, -2, 2), axis=2),
        "mean": s.mean(axis=2), "median": np.median(s, axis=2),
        "max": s.max(axis=2), "unbounded": logsumexp(s, axis=2),
    }
    result = {}
    for name, p in pools.items():
        d = p[:, 1] - p[:, 0]
        result[name] = float(np.mean((d > 0) + .5 * (d == 0)))

    def med(x):
        m = np.median(x, axis=2)
        return m[:, 1] - m[:, 0]

    result["median_expression"] = concordance(med(a), med(b), "block_logistic")
    result["pooled_counts"] = concordance(train_features["pooled_logistic"],
                                          test_features["pooled_logistic"], "pooled_logistic")
    return result


def main():
    rng = np.random.default_rng(20261123)
    regimes = [("localized_reversal", False), ("localized_reversal", True),
               ("diffuse_reversal", True), ("null", True)]
    results, checks = {}, {}
    for regime, contamination in regimes:
        rows = []
        for _ in range(100):
            train = draw(rng, 24, "localized_reversal")
            test = draw(rng, 80, regime)
            if contamination:
                corrupt(rng, test)
            rows.append(run(train, test))
        means = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
        key = regime + ("_corrupt" if contamination else "_clean")
        results[key] = means
        gap = means["bounded"] - max(v for k, v in means.items() if k != "bounded")
        if regime == "null":
            checks[key] = abs(means["bounded"] - .5) <= .02
        else:
            checks[key] = gap >= (.01 if key == "localized_reversal_corrupt" else -.02)
    path = Path(__file__).resolve().parents[1] / "output/trackA_pooling_simulation/result.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"label": "simulation", "passed": all(checks.values()),
                                "checks": checks, "concordance": results}, indent=2) + "\n",
                    encoding="utf-8")


if __name__ == "__main__":
    main()
