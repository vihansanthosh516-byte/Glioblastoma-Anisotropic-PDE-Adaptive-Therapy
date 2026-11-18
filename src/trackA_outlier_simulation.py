import json
from pathlib import Path

import numpy as np
from trackA_order_simulation import draw, contrasts, concordance


def features(rng, n, regime, contaminated):
    x = draw(rng, n, regime)
    if contaminated:
        patients = np.arange(n)[:, None]
        regions = np.arange(2)
        blocks = rng.integers(4, size=(n, 2))
        original = x[patients, regions, blocks]
        adjusted = original * np.exp(rng.normal(scale=3, size=original.shape))
        probability = adjusted / adjusted.sum(axis=-1, keepdims=True)
        x[patients, regions, blocks] = rng.multinomial(original.sum(axis=-1), probability)
    arms = contrasts(x)
    median = np.median(np.log1p(10000 * x / x.sum(axis=-1, keepdims=True)), axis=2)
    arms["median_logistic"] = median[:, 1] - median[:, 0]
    arms["combined_logistic"] = np.column_stack((arms["median_logistic"], arms["soft"]))
    return arms


def main():
    rng = np.random.default_rng(20261121)
    results = {}
    candidate = "combined_logistic"
    for regime in ("diffuse_reversal", "null"):
        for contaminated in (False, True):
            runs = []
            for _ in range(100):
                train = features(rng, 24, regime, contaminated)
                test = features(rng, 80, regime, contaminated)
                runs.append([concordance(train[k], test[k], k) for k in train])
            suffix = "contaminated" if contaminated else "clean"
            results[f"{regime}_{suffix}"] = dict(zip(train, np.mean(runs, axis=0).tolist()))
    gaps = {k: d[candidate] - max(v for m, v in d.items() if m != candidate)
            for k, d in results.items()}
    checks = {
        "contaminated_superiority": gaps["diffuse_reversal_contaminated"] >= .01,
        "clean_tolerance": gaps["diffuse_reversal_clean"] >= -.02,
        **{k: abs(d[candidate] - .50) <= .02 for k, d in results.items() if k.startswith("null")},
    }
    path = Path(__file__).resolve().parents[1] / "output/trackA_outlier_simulation/result.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"label": "simulation", "passed": all(checks.values()),
                                "checks": checks, "concordance": results, "gaps": gaps}, indent=2) + "\n",
                    encoding="utf-8")


if __name__ == "__main__":
    main()
