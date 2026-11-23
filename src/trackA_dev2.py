"""Track A development falsification 2 on generated data only.

Plan: paper/trackA_dev2_falsification_plan_2026-11-23.md (frozen before run).
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from .trackA_regional_model import NONSPATIAL, fit_regional
from .trackA_tmz_model import Cultures, _fit_shared, _groups, _means, _shuffle_links

DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "output/trackA_dev2/result.json"
SEED = 202611232
REPLICATES = 20
GRID = tuple(np.round(np.arange(0, 3.01, 0.25), 2))
MIN_OVERRIDES, HARM_CEILING, LEVEL = 5, 0.35, 0.95
PENALTY = 0.2


# ---------------------------------------------------------------- regional
def regional_fixture(rng, count, *, density_only=False):
    orientation = rng.choice([-1, 1], size=count)
    nonspatial = rng.normal(size=(count, 3)) + [0.35, 0.15, 0.45]
    if density_only:
        spatial = nonspatial[:, 1] + rng.normal(scale=0.3, size=count)
    else:
        spatial = rng.normal(loc=0.9, scale=0.7, size=count)
    layers = {
        "composition": nonspatial[:, 0], "density_depth": nonspatial[:, 1],
        "malignant_expression": nonspatial[:, 2],
        "nominal_spatial": spatial + rng.normal(scale=0.08, size=count),
        "alternative_spatial": spatial + rng.normal(scale=0.08, size=count),
        "molecule_spatial": spatial + rng.normal(scale=1.2, size=count),
    }
    return {k: (v * orientation)[:, None] for k, v in layers.items()}, orientation


def offset(features, orientation, nominal, alternative, mask=None):
    mask = np.ones(len(orientation), bool) if mask is None else mask
    shift = orientation * mask
    return {**features,
            "nominal_spatial": features["nominal_spatial"] + nominal * shift[:, None],
            "alternative_spatial": features["alternative_spatial"] + alternative * shift[:, None]}


def proposals(model, features):
    scores, _ = model.predict(features, bound=np.inf,
                              audit_qualified=np.zeros(len(features["composition"]), bool))
    base = scores["nonspatial"]
    rec = model.reconstruction_margins(features)
    agree = (np.sign(rec[:, 0]) == np.sign(rec[:, 1])) & (rec[:, 0] != 0)
    proposed = agree & (np.sign(rec[:, 0]) != np.sign(base))
    return scores, base, rec, proposed, np.min(np.abs(rec), axis=1)


def cp_upper(wrong, n):
    return 1.0 if wrong >= n else float(beta.ppf(LEVEL, wrong + 1, n - wrong))


def harm_ok(proposed, strength, wrong, threshold):
    selected = proposed & (strength >= threshold)
    n = int(selected.sum())
    return n >= MIN_OVERRIDES and cp_upper(int((selected & wrong).sum()), n) <= HARM_CEILING


def choose_threshold(model, split_a, split_b):
    """Choose on split A, validate on split B. Return None when rejected."""
    rule = {}
    for name, (features, orientation) in (("a", split_a), ("b", split_b)):
        _, _, rec, proposed, strength = proposals(model, features)
        rule[name] = (proposed, strength, rec[:, 0] * orientation <= 0)
    chosen = next((t for t in GRID if harm_ok(*rule["a"], t)), None)
    if chosen is None or not harm_ok(*rule["b"], chosen):
        return None
    return float(chosen)


def correct(margin, orientation):
    product = margin * orientation
    return np.where(product > 0, 1.0, np.where(product == 0, 0.5, 0.0))


def regional_replicate(rng, scenario):
    density_only = scenario == "density_only"
    draw = lambda: regional_fixture(rng, 160, density_only=density_only)
    train, cal_a, cal_b, ev = draw(), draw(), draw(), draw()
    model = fit_regional(train[0], train[1], strength=0.1)
    shift = {"clean": (0, 0), "density_only": (0, 0), "shared_error": (-1.8, -1.8),
             "divergent_error": (-1.8, 1.8), "calibration_mismatch": (-1.8, -1.8),
             "missed_contamination": (-1.8, -1.8)}[scenario]
    clean_cal = scenario in ("calibration_mismatch", "missed_contamination")
    cals = [c if clean_cal else (offset(c[0], c[1], *shift), c[1]) for c in (cal_a, cal_b)]
    mask = rng.random(160) < 0.3 if scenario == "missed_contamination" else None
    measured = offset(ev[0], ev[1], *shift, mask)
    threshold = choose_threshold(model, *cals)
    scores, base, rec, proposed, strength = proposals(model, measured)
    used = proposed & (strength >= threshold) if threshold is not None else np.zeros_like(proposed)
    scores["override_calibrated"] = np.where(used, rec[:, 0], base)
    scores["override_unchecked"] = np.where(proposed, rec[:, 0], base)
    hit = {arm: correct(s, ev[1]) for arm, s in scores.items()}
    harmed = used & (hit["override_calibrated"] < hit["nonspatial"])
    improved = used & (hit["override_calibrated"] > hit["nonspatial"])
    unchecked_harm = proposed & (hit["override_unchecked"] < hit["nonspatial"])
    return {
        "threshold": threshold,
        "concordance": {arm: float(h.mean()) for arm, h in hit.items()},
        "proposed": int(proposed.sum()), "overrides": int(used.sum()),
        "harmed": int(harmed.sum()), "improved": int(improved.sum()),
        "unchecked_harmed": int(unchecked_harm.sum()),
        "per_patient": [
            {"baseline": float(b), "nominal": float(r[0]), "alternative": float(r[1]),
             "final": float(f), "changed": bool(u), "harmed": bool(h)}
            for b, r, f, u, h in zip(base, rec, scores["override_calibrated"], used, harmed, strict=True)
        ],
    }


# --------------------------------------------------------------------- TMZ
W = np.array([0.9, -0.6, 0, 0, 0, 0, 0, 0])


def tmz_cultures(rng, prefix, count, *, between, within, max_size=4):
    sizes = rng.integers(1, max_size + 1, size=count)
    groups = np.repeat(np.arange(count), sizes)
    common, noise = rng.normal(size=(count, 8)), rng.normal(size=(len(groups), 8))
    x = 0.6 * common[groups] + noise
    y = (0.6 * common @ between)[groups] + noise @ within
    y += rng.normal(scale=0.3, size=count)[groups] + rng.normal(scale=0.3, size=len(groups))
    return Cultures(x, y, [f"{prefix}-{g:03d}" for g in groups],
                    [f"{prefix}-c{i:04d}" for i in range(len(groups))])


def within_direction(x, y, origins):
    _, inv, sizes = _groups(origins)
    keep = sizes[inv] > 1
    if not keep.any():
        return None
    cx = x - _means(x, inv, sizes)[inv]
    cy = y - _means(y, inv, sizes)[inv]
    d = Ridge(alpha=PENALTY, fit_intercept=False).fit(cx[keep], cy[keep]).coef_
    return d if np.any(d) else None


def blind_direction(x, origins):
    _, inv, sizes = _groups(origins)
    keep = sizes[inv] > 1
    cx = (x - _means(x, inv, sizes)[inv])[keep]
    return np.linalg.svd(cx, full_matrices=False)[2][0]


def balanced_fit(x, origins, y):
    _, inv, counts = _groups(origins)
    return Ridge(alpha=PENALTY).fit(x, y, sample_weight=1.0 / (len(counts) * counts[inv]))


def origin_mae(pred, y, origins):
    _, inv, counts = _groups(origins)
    return float(np.abs(_means(pred, inv, counts) - _means(y, inv, counts)).mean())


def augmented(x, origins, y, d, eval_x):
    model = balanced_fit(np.column_stack((x, x @ d)), origins, y)
    return model.predict(np.column_stack((eval_x, eval_x @ d)))


def admit_direction(x, y, origins, rng):
    names = rng.permutation(np.unique(origins))
    half = np.isin(origins, names[: len(names) // 2])
    sides = (half, ~half)
    dirs = [within_direction(x[s], y[s], origins[s]) for s in sides]
    if any(d is None for d in dirs):
        return False, None
    cosine = float(dirs[0] @ dirs[1] / np.linalg.norm(dirs[0]) / np.linalg.norm(dirs[1]))
    if cosine <= 0:
        return False, cosine
    for fit_side, test_side, d in ((sides[0], sides[1], dirs[0]), (sides[1], sides[0], dirs[1])):
        base = balanced_fit(x[fit_side], origins[fit_side], y[fit_side]).predict(x[test_side])
        aug = augmented(x[fit_side], origins[fit_side], y[fit_side], d, x[test_side])
        if origin_mae(aug, y[test_side], origins[test_side]) >= origin_mae(base, y[test_side], origins[test_side]):
            return False, cosine
    return True, cosine


def tmz_replicate(rng, scenario):
    between = {"transfer": W, "reversal": -W, "unequal_sisters": W,
               "assay_shift": W, "null": 0 * W}[scenario]
    within = 0 * W if scenario == "null" else W
    size = 8 if scenario == "unequal_sisters" else 4
    train = tmz_cultures(rng, "t", 40, between=between, within=within, max_size=size)
    ev = tmz_cultures(rng, "e", 100, between=between, within=within, max_size=size)
    eval_raw = ev.expression.copy()
    if scenario == "assay_shift":
        eval_raw[:, 2] += 1.5 * eval_raw[:, 0]
    _, inv, counts = _groups(train.origins)
    scale = StandardScaler().fit(train.expression, sample_weight=1.0 / (len(counts) * counts[inv]))
    x, ex, y, o = scale.transform(train.expression), scale.transform(eval_raw), train.viability, train.origins
    shuffled = _shuffle_links(o, ["assay"] * len(o), int(rng.integers(1 << 31)))
    preds = {"origin_balanced_ridge": balanced_fit(x, o, y).predict(ex)}
    admitted, cosine = admit_direction(x, y, o, np.random.default_rng(int(rng.integers(1 << 31))))
    d_true = within_direction(x, y, o)
    preds["replicated_direction"] = augmented(x, o, y, d_true, ex) if admitted else preds["origin_balanced_ridge"]
    preds["direction_unchecked"] = augmented(x, o, y, d_true, ex)
    preds["shuffled_link_direction"] = augmented(x, o, y, within_direction(x, y, shuffled), ex)
    preds["outcome_blind_direction"] = augmented(x, o, y, blind_direction(x, o), ex)
    intercept, slope = _fit_shared(x, y, o, o, PENALTY, 1.0)
    preds["one_slope_contrast"] = intercept + ex @ slope
    mae = {arm: origin_mae(p, ev.viability, ev.origins) for arm, p in preds.items()}
    return {"admitted": bool(admitted), "cosine": cosine, "origin_mae": mae}


# ------------------------------------------------------------------ driver
def summarize(values):
    values = np.asarray(values, float)
    half = 1.96 * values.std(ddof=1) / np.sqrt(len(values))
    return {"mean": float(values.mean()), "ci95": [float(values.mean() - half), float(values.mean() + half)]}


def run():
    root = np.random.SeedSequence(SEED)
    reg_seq, tmz_seq = root.spawn(2)
    regional = {}
    reg_scen = ("clean", "shared_error", "divergent_error", "calibration_mismatch",
                "density_only", "missed_contamination")
    for scen, seq in zip(reg_scen, reg_seq.spawn(len(reg_scen)), strict=True):
        reps = [regional_replicate(np.random.default_rng(s), scen) for s in seq.spawn(REPLICATES)]
        arms = reps[0]["concordance"].keys()
        others = [a for a in arms if a not in ("override_calibrated", "override_unchecked")]
        strongest = max(others, key=lambda a: np.mean([r["concordance"][a] for r in reps]))
        regional[scen] = {
            "mean_concordance": {a: float(np.mean([r["concordance"][a] for r in reps])) for a in arms},
            "strongest_noncandidate": strongest,
            "candidate_minus_strongest": summarize(
                [r["concordance"]["override_calibrated"] - r["concordance"][strongest] for r in reps]),
            "candidate_minus_unchecked": summarize(
                [r["concordance"]["override_calibrated"] - r["concordance"]["override_unchecked"] for r in reps]),
            "harms_per_patient": float(np.mean([r["harmed"] / 160 for r in reps])),
            "harms_per_override": float(sum(r["harmed"] for r in reps) / max(1, sum(r["overrides"] for r in reps))),
            "improved_total": int(sum(r["improved"] for r in reps)),
            "harmed_total": int(sum(r["harmed"] for r in reps)),
            "overrides_total": int(sum(r["overrides"] for r in reps)),
            "unchecked_harmed_total": int(sum(r["unchecked_harmed"] for r in reps)),
            "thresholds": [r["threshold"] for r in reps],
            "replicates": reps,
        }
    tmz = {}
    tmz_scen = ("transfer", "reversal", "unequal_sisters", "assay_shift", "null")
    for scen, seq in zip(tmz_scen, tmz_seq.spawn(len(tmz_scen)), strict=True):
        reps = [tmz_replicate(np.random.default_rng(s), scen) for s in seq.spawn(REPLICATES)]
        arms = reps[0]["origin_mae"].keys()
        tmz[scen] = {
            "mean_origin_mae": {a: float(np.mean([r["origin_mae"][a] for r in reps])) for a in arms},
            "base_minus_candidate": summarize([r["origin_mae"]["origin_balanced_ridge"]
                                               - r["origin_mae"]["replicated_direction"] for r in reps]),
            "shuffled_minus_candidate": summarize([r["origin_mae"]["shuffled_link_direction"]
                                                   - r["origin_mae"]["replicated_direction"] for r in reps]),
            "admitted_fraction": float(np.mean([r["admitted"] for r in reps])),
            "replicates": reps,
        }
    adverse = ("shared_error", "divergent_error", "calibration_mismatch", "density_only", "missed_contamination")
    rules = {
        "regional_clean_beats_strongest": regional["clean"]["candidate_minus_strongest"]["ci95"][0] > 0,
        "regional_harm_controlled": all(regional[s]["harms_per_patient"] <= 0.02 for s in adverse),
        "regional_safeguard_identifiable": any(regional[s]["candidate_minus_unchecked"]["ci95"][0] > 0 for s in adverse),
        "tmz_transfer_beats_base": tmz["transfer"]["base_minus_candidate"]["ci95"][0] > 0,
        "tmz_transfer_beats_shuffle": tmz["transfer"]["shuffled_minus_candidate"]["ci95"][0] > 0,
        "tmz_reversal_safe": tmz["reversal"]["base_minus_candidate"]["mean"] >= -0.01,
        "tmz_null_safe": tmz["null"]["base_minus_candidate"]["mean"] >= -0.01,
    }
    return {"scope": "generated development data only", "plan": "paper/trackA_dev2_falsification_plan_2026-11-23.md",
            "seed": SEED, "replicates": REPLICATES, "rules": rules, "regional": regional, "tmz": tmz}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_suffix(".partial")
    with tmp.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=1, allow_nan=False)
        stream.write("\n")
    tmp.rename(args.output)  # fails if output already exists on Windows: no silent rerun
    print(json.dumps(result["rules"], indent=1))


if __name__ == "__main__":
    main()
