"""Frozen Track A synthetic spatial stress simulation."""
import json
from pathlib import Path


def simulate():
    state = 20261124

    def random():
        nonlocal state
        state = (1664525 * state + 1013904223) % 2**32
        return state / 2**32

    def bridge_score(labels, edges, positions):
        if not edges:
            return 0.0, 0.0
        observed = sum(labels[i] != labels[j] for i, j in edges) / len(edges)
        n = len(positions)
        a = sum(labels[i] for i in positions)
        return observed, observed - 2 * a * (n - a) / (n * (n - 1))

    results = {}
    for regime in ("clean_signal", "noisy_signal", "noisy_null"):
        totals = dict.fromkeys(("depth_matched", "global_adjusted", "raw", "composition", "degree_composition", "geometry"), 0.0)
        empty = [0, 0]
        for _ in range(1000):
            positions = list(range(24))
            graph = [(i, (i + 1) % 24) for i in positions]
            graph += [(i, (i + 3) % 24) for i in positions if random() < .25]
            offset = int(random() * 24)
            regions = []
            for region, depth_probability in enumerate((.25, .75)):
                low = [random() < depth_probability for _ in positions]
                latent = [int((i + offset) % 24 < 12) if region == 0 or regime == "noisy_null"
                          else (i + offset) % 2 for i in positions]
                labels = [v ^ (random() < (.20 if shallow else .02)) if regime != "clean_signal" else v
                          for v, shallow in zip(latent, low)]
                high = [i for i in positions if not low[i]]
                eligible = [(i, j) for i, j in graph if not low[i] and not low[j]]
                empty[region] += not eligible
                raw, adjusted = bridge_score(labels, graph, positions)
                _, matched = bridge_score(labels, eligible, high)
                weighted = sum(labels[i] + labels[j] for i, j in graph) / (2 * len(graph))
                regions.append((matched, adjusted, raw, sum(labels) / 24, weighted, 0))
            for key, core, edge in zip(totals, *regions):
                totals[key] += 1 if edge > core else .5 if edge == core else 0
        results[regime] = {"concordance": {k: v / 1000 for k, v in totals.items()},
                           "no_high_depth_bridge_rate": dict(zip(("core", "edge"), (n / 1000 for n in empty)))}
    noisy = results["noisy_signal"]["concordance"]
    checks = {"noisy_superiority": noisy["depth_matched"] - max(v for k, v in noisy.items() if k != "depth_matched") >= .01,
              "clean": results["clean_signal"]["concordance"]["depth_matched"] >= .80,
              "null": abs(results["noisy_null"]["concordance"]["depth_matched"] - .50) <= .03}
    return {"label": "simulation", "passed": all(checks.values()), "checks": checks, "regimes": results}


if __name__ == "__main__":
    path = Path(__file__).resolve().parents[1] / "output/trackA_spatial_stress/result.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(simulate(), indent=2) + "\n", encoding="utf-8")
