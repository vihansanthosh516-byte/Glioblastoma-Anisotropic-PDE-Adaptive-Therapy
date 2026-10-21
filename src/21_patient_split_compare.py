"""
21_patient_split_compare.py
===========================
Paired patient-level comparison for the donor-disjoint test set (5 patients).

Decision rule (fixed in research_notes/plan_patient_split_methods_2026-10-19.md):
  A new method "helps" only if
    (a) its seed-mean test macro-F1 exceeds the best baseline's seed-mean test macro-F1, AND
    (b) the 95% bootstrap interval of (method - best baseline) excludes zero.
  Otherwise the result is "no measurable gain".

Bootstrap: resample the 5 test patients with replacement (2000 times). Each resample
takes all cells of the chosen patients; macro-F1 is computed for both models on the same
cells, and the difference is recorded. Seed 42 predictions are used for the interval;
seed means are used for the point estimate.

Inputs:
  output/patient_split_seeds/seed{S}/test_predictions_*_seed{S}.npz   (baselines, from the seed sweep)
  output/patient_split_v2_methods/pred_seed{S}_*.npz                  (new methods, from script 20)
Output:
  output/patient_split_v2_methods/comparison.json

Run: python src/21_patient_split_compare.py
"""
import os, glob, json
import numpy as np
from sklearn.metrics import f1_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "output")
SEEDS = [42, 43, 44, 45, 46]
PID = np.load(os.path.join(OUT, "nn_patient_ids.npy"), allow_pickle=True)
N_BOOT = 2000

BASELINES = {
    "attention_net_05": "patient_split_seeds/seed{S}/test_predictions_nn_seed{S}.npz",
    "transformer_07": "patient_split_seeds/seed{S}/test_predictions_transformer_seed{S}.npz",
    "hybrid_08": "patient_split_seeds/seed{S}/test_predictions_hybrid_seed{S}.npz",
    "gat_13": "patient_split_seeds/seed{S}/test_predictions_gat13_seed{S}.npz",
    "cgat_16": "patient_split_seeds/seed{S}/test_predictions_cgat16_seed{S}.npz",
}
NEW = {
    "M1p_pca_lr": "patient_split_v2_methods/pred_seed{S}_M1p_pca_lr.npz",
    "M1p_pca_rf": "patient_split_v2_methods/pred_seed{S}_M1p_pca_rf.npz",
    "M3p_panel_lr": "patient_split_v2_methods/pred_seed{S}_M3p_panel_lr.npz",
}


def load(rel, S):
    p = os.path.join(OUT, rel.format(S=S))
    if not os.path.exists(p):
        return None
    d = np.load(p)
    return d["y_true"], d["y_pred"], d["test_index"]


def seed_table(spec):
    """Return {seed: (y_true, y_pred, test_index)} for seeds where the file exists."""
    return {S: r for S in SEEDS if (r := load(spec, S)) is not None}


def macro(yt, yp):
    return f1_score(yt, yp, average="macro")


def bootstrap_diff(a, b, rng):
    """a, b: (y_true, y_pred, test_index) on the same cells. Patient-level bootstrap of macro-F1(a) - macro-F1(b)."""
    yt = a[0]
    pats = PID[a[2]]
    uniq = np.unique(pats)
    diffs = []
    for _ in range(N_BOOT):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([np.where(pats == p)[0] for p in pick])
        diffs.append(macro(yt[idx], a[1][idx]) - macro(b[0][idx], b[1][idx]))
    return np.percentile(diffs, [2.5, 97.5]).tolist()


def main():
    base = {n: seed_table(s) for n, s in BASELINES.items()}
    base = {n: t for n, t in base.items() if t}
    new = {n: seed_table(s) for n, s in NEW.items()}
    new = {n: t for n, t in new.items() if t}

    # seed-mean macro-F1 per model
    mean_f1 = {}
    for name, tab in {**base, **new}.items():
        vals = [macro(yt, yp) for (yt, yp, _) in tab.values()]
        mean_f1[name] = {"seed_mean": float(np.mean(vals)), "min": float(np.min(vals)),
                         "max": float(np.max(vals)), "n_seeds": len(vals)}

    best_name = max(base, key=lambda n: mean_f1[n]["seed_mean"])
    best_base = base[best_name]
    rng = np.random.default_rng(42)
    results = {"best_baseline": best_name, "mean_macro_f1": mean_f1, "comparisons": {}}
    for name, tab in new.items():
        if 42 not in tab or 42 not in best_base:
            results["comparisons"][name] = {"status": "missing seed-42 predictions"}
            continue
        ci = bootstrap_diff(tab[42], best_base[42], rng)
        diff_mean = mean_f1[name]["seed_mean"] - mean_f1[best_name]["seed_mean"]
        helps = (diff_mean > 0) and (ci[0] > 0)
        results["comparisons"][name] = {
            "diff_seed_mean": diff_mean, "bootstrap_95ci_seed42": ci,
            "verdict": "helps" if helps else "no measurable gain"}

    dest = os.path.join(OUT, "patient_split_v2_methods", "comparison.json")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "w") as f:
        json.dump(results, f, indent=2)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
