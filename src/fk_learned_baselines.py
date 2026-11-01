"""Fold-trained U-Net mask baseline and training-only Poisson count model for Track B.

Both fit on outer-training patients only, inside the outer-fold nuisance fit. Inputs are the visit-1 core, the
visit-2 core and the origin-only domain; the visit-3 core is a training target only. The count model forecasts a
rate per day from history; the scorer multiplies it by the actual interval, so the forecast never sees the target
day. The count model is never turned into a mask and never scored with Dice. Settings are frozen in
``paper/track_b_fixed_growth_baseline_spec_2026-10-31.md`` section 5.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, asdict
import hashlib

import numpy as np
import torch
from scipy.ndimage import distance_transform_edt, label

from src.fk_mu_scoring import dice_with_empty_rule


@dataclass(frozen=True)
class UNetRecipe:
    widths: tuple[int, ...] = (16, 32, 64, 128)
    learning_rate: float = 1e-3
    betas: tuple[float, float] = (0.9, 0.999)
    weight_decay: float = 0.0
    epochs: int = 300
    n_seeds: int = 3
    seed_offsets: tuple[int, ...] = (0, 1000, 2000)
    base_seed: int = 20261102
    threshold: float = 0.5
    dice_smoothing: float = 1.0
    leaky_slope: float = 0.01
    # None = frozen full-grid recipe. An integer is the PROPOSED, NOT ADOPTED history-centred cube crop.
    crop_voxels: int | None = None


UNET_RECIPE = UNetRecipe()
NEW_PIECE_MARGIN_MM = 2.0
NEW_PIECE_MIN_VOXELS = 2
COUNT_MAX_ITER = 50
COUNT_TOL = 1e-8


def _mask(value, *, name: str) -> np.ndarray:
    array = np.asarray(value)
    if array.ndim != 3 or (array.dtype != np.bool_ and not np.isin(array, (0, 1)).all()):
        raise ValueError(f"{name} must be a 3D binary mask")
    return array.astype(bool, copy=False)


# ---------------------------------------------------------------- U-Net

def _block(c_in, c_out, slope):
    return torch.nn.Sequential(
        torch.nn.Conv3d(c_in, c_out, 3, padding=1), torch.nn.InstanceNorm3d(c_out, affine=True),
        torch.nn.LeakyReLU(slope),
        torch.nn.Conv3d(c_out, c_out, 3, padding=1), torch.nn.InstanceNorm3d(c_out, affine=True),
        torch.nn.LeakyReLU(slope),
    )


class UNet3D(torch.nn.Module):
    def __init__(self, widths=UNET_RECIPE.widths, slope=UNET_RECIPE.leaky_slope, in_channels=3):
        super().__init__()
        self.down = torch.nn.ModuleList()
        c = in_channels
        for w in widths:
            self.down.append(_block(c, w, slope))
            c = w
        self.pool = torch.nn.MaxPool3d(2)
        self.up = torch.nn.ModuleList()
        self.up_blocks = torch.nn.ModuleList()
        for w in reversed(widths[:-1]):
            self.up.append(torch.nn.ConvTranspose3d(c, w, 2, stride=2))
            self.up_blocks.append(_block(2 * w, w, slope))
            c = w
        self.head = torch.nn.Conv3d(c, 1, 1)

    def forward(self, x):
        skips = []
        for i, block in enumerate(self.down):
            x = block(x)
            if i < len(self.down) - 1:
                skips.append(x)
                x = self.pool(x)
        for up, block, skip in zip(self.up, self.up_blocks, reversed(skips)):
            x = block(torch.cat([up(x), skip], dim=1))
        return self.head(x)


def _padded_shape(shape, levels):
    m = 2 ** (levels - 1)
    padded = [int(np.ceil(s / m) * m) for s in shape]
    padded[0] = padded[1] = max(padded[0], padded[1])
    return tuple(padded)


def _pad(array, shape):
    out = np.zeros(shape, dtype=np.float32)
    out[tuple(slice(0, s) for s in array.shape)] = array
    return out


def unet_inputs(visit1, visit2, domain) -> np.ndarray:
    """Stack the three binary history channels. No target-derived data enters."""
    return np.stack([_mask(visit1, name="visit1").astype(np.float32), _mask(visit2, name="visit2").astype(np.float32),
                     _mask(domain, name="domain").astype(np.float32)])


def _loss(logits, target, smoothing):
    prob = torch.sigmoid(logits)
    inter = (prob * target).sum()
    dice = 1 - (2 * inter + smoothing) / (prob.sum() + target.sum() + smoothing)
    return dice + torch.nn.functional.binary_cross_entropy_with_logits(logits, target)


def train_unet_models(pairs, *, fold: int, recipe: UNetRecipe = UNET_RECIPE) -> dict:
    """Train ``recipe.n_seeds`` models on (inputs[3,X,Y,Z], target[X,Y,Z]) pairs from outer-training patients."""
    if not pairs:
        return {"status": "unavailable", "reason": "no usable training pairs", "models": []}
    shape = pairs[0][1].shape
    padded = _padded_shape(shape, len(recipe.widths))
    xs = [torch.from_numpy(np.stack([_pad(c, padded) for c in x]))[None] for x, _ in pairs]
    ys = [torch.from_numpy(_pad(y.astype(np.float32), padded))[None, None] for _, y in pairs]
    prior_det = torch.are_deterministic_algorithms_enabled()
    torch.use_deterministic_algorithms(True)
    models, seeds = [], []
    try:
        for offset in recipe.seed_offsets[:recipe.n_seeds]:
            seed = recipe.base_seed + int(fold) + int(offset)
            seeds.append(seed)
            torch.manual_seed(seed)
            gen = np.random.default_rng(seed)
            model = UNet3D(recipe.widths, recipe.leaky_slope)
            opt = torch.optim.Adam(model.parameters(), lr=recipe.learning_rate, betas=recipe.betas,
                                   weight_decay=recipe.weight_decay)
            model.train()
            for _ in range(recipe.epochs):
                for i in gen.permutation(len(xs)):
                    x, y = xs[i], ys[i]
                    for axis in (2, 3, 4):
                        if gen.random() < 0.5:
                            x, y = torch.flip(x, (axis,)), torch.flip(y, (axis,))
                    k = int(gen.integers(4))
                    x, y = torch.rot90(x, k, (2, 3)), torch.rot90(y, k, (2, 3))
                    opt.zero_grad()
                    loss = _loss(model(x), y, recipe.dice_smoothing)
                    if not torch.isfinite(loss):
                        return {"status": "unavailable", "reason": "non-finite loss", "models": [], "seeds": seeds}
                    loss.backward()
                    opt.step()
            model.eval()
            models.append(model)
    finally:
        torch.use_deterministic_algorithms(prior_det)
    result = {"status": "ok", "models": models, "seeds": seeds, "padded_shape": padded, "shape": shape,
              "recipe": asdict(recipe), "n_training_pairs": len(pairs)}
    if not any(predict_unet(result, x, domain=np.ones(shape, bool))["mask"].any()
               for x, _ in pairs):
        return {"status": "unavailable", "reason": "all training predictions empty", "models": [], "seeds": seeds}
    return result


def predict_unet(trained: Mapping, inputs, *, domain) -> dict:
    """Average the seed probability maps, threshold at the frozen value, clip to the origin-only domain."""
    if trained.get("status") != "ok":
        raise ValueError("U-Net is unavailable for this fold")
    shape = tuple(trained["shape"])
    x = torch.from_numpy(np.stack([_pad(c, trained["padded_shape"]) for c in inputs]))[None]
    with torch.no_grad():
        probs = [torch.sigmoid(m(x))[0, 0].numpy()[tuple(slice(0, s) for s in shape)] for m in trained["models"]]
    mean = np.mean(probs, axis=0)
    mask = (mean >= trained["recipe"]["threshold"]) & _mask(domain, name="domain")
    return {"mask": mask, "probability": mean, "seed_probabilities": probs}


def history_crop_box(history_core, shape, size: int) -> tuple[slice, ...]:
    """Cube of ``size`` voxels centred on the history core's centroid, shifted to lie inside the grid.

    Only the visit-2 (history) core sets the box. Axes shorter than ``size`` use the whole axis.
    """
    core = _mask(history_core, name="history core")
    if not core.any():
        raise ValueError("crop needs a nonempty history core")
    centre = np.rint(np.argwhere(core).mean(axis=0)).astype(int)
    box = []
    for c, n in zip(centre, shape):
        width = min(int(size), int(n))
        start = int(np.clip(c - width // 2, 0, n - width))
        box.append(slice(start, start + width))
    return tuple(box)


def uncrop_mask(mask, box, shape) -> np.ndarray:
    full = np.zeros(shape, dtype=bool)
    full[box] = mask
    return full


def unet_training_fingerprint(patient_ids, pairs, *, fold: int, recipe: UNetRecipe) -> str:
    """SHA-256 over everything that determines the trained U-Net models (inputs, targets, fold, recipe)."""
    h = hashlib.sha256()
    h.update(repr((int(fold), sorted(asdict(recipe).items()))).encode())
    for pid, (x, y) in zip(patient_ids, pairs):
        h.update(str(pid).encode())
        h.update(np.ascontiguousarray(x, dtype=np.float32).tobytes())
        h.update(np.packbits(np.asarray(y, bool)).tobytes())
        h.update(repr((tuple(x.shape), tuple(y.shape))).encode())
    return h.hexdigest()


# ---------------------------------------------------------------- count model

def new_label_pieces(previous, current, *, spacing_mm, domain=None) -> int:
    """Count 26-connected pieces of current core outside previous core grown by 2 mm (>= 2 voxels)."""
    prev = _mask(previous, name="previous core")
    cur = _mask(current, name="current core")
    if domain is not None:
        d = _mask(domain, name="domain")
        prev, cur = prev & d, cur & d
    grown = (distance_transform_edt(~prev, sampling=spacing_mm) <= NEW_PIECE_MARGIN_MM + 1e-9) if prev.any() \
        else np.zeros_like(prev)
    labels, n = label(cur & ~grown, structure=np.ones((3, 3, 3), bool))
    if n == 0:
        return 0
    sizes = np.bincount(labels.ravel())[1:]
    return int(np.count_nonzero(sizes >= NEW_PIECE_MIN_VOXELS))


def fit_count_model(counts, log_volumes, intervals_days) -> dict:
    """Poisson GLM log mu = b0 + b1 log V + log t by Newton-IRLS; mean-rate fallback on failure."""
    y = np.asarray(counts, float)
    t = np.asarray(intervals_days, float)
    if y.size == 0 or np.any(t <= 0) or np.any(y < 0):
        raise ValueError("count model needs nonnegative counts and positive intervals")
    # A zero training total would give rate 0 and infinite deviance for any new piece; add 0.5 to the count.
    mean_rate = float((y.sum() if y.sum() > 0 else 0.5) / t.sum())
    out = {"mean_rate_per_day": mean_rate, "n_training": int(y.size)}
    lv = np.asarray(log_volumes, float)
    usable = np.isfinite(lv)
    if usable.sum() < 3 or np.ptp(lv[usable]) == 0:
        return {**out, "status": "mean_rate_fallback", "reason": "fewer than 3 usable patients or constant volume"}
    X = np.column_stack([np.ones(usable.sum()), lv[usable]])
    yy, off = y[usable], np.log(t[usable])
    beta = np.array([np.log(max(yy.sum(), 0.5) / t[usable].sum()), 0.0])
    for it in range(COUNT_MAX_ITER):
        mu = np.exp(X @ beta + off)
        H = X.T @ (mu[:, None] * X)
        try:
            step = np.linalg.solve(H, X.T @ (yy - mu))
        except np.linalg.LinAlgError:
            return {**out, "status": "mean_rate_fallback", "reason": "singular step"}
        beta = beta + step
        if not np.isfinite(beta).all():
            return {**out, "status": "mean_rate_fallback", "reason": "non-finite coefficients"}
        if np.max(np.abs(step)) < COUNT_TOL:
            return {**out, "status": "glm", "coef": beta.tolist(), "iterations": it + 1}
    return {**out, "status": "mean_rate_fallback", "reason": "no convergence"}


def count_rate_per_day(model: Mapping, log_volume: float) -> float:
    if model["status"] == "glm" and np.isfinite(log_volume):
        b0, b1 = model["coef"]
        return float(np.exp(b0 + b1 * log_volume))
    return float(model["mean_rate_per_day"])


def count_metrics(mu: float, y: int) -> dict:
    mu = float(mu)
    if not np.isfinite(mu) or mu <= 0:
        raise ValueError("predicted count must be positive and finite")
    deviance = 2 * ((y * np.log(y / mu) if y > 0 else 0.0) - (y - mu))
    return {"poisson_deviance": float(deviance), "abs_error": float(abs(y - mu)),
            "brier_any_new": float(((1 - np.exp(-mu)) - float(y >= 1)) ** 2)}


# ---------------------------------------------------------------- runner wiring

def add_learned_baselines_arm(fit_nuisance, forecast, fallback_forecast, score, *, domain_by_patient: Mapping,
                              spacing_mm, valid_support_by_patient: Mapping | None = None,
                              recipe: UNetRecipe = UNET_RECIPE):
    """Wrap patient-CV callbacks to add the fold-trained U-Net mask arm and the count model.

    ``fit_nuisance`` is called with outer-training records only. Training uses visits 1-2 as inputs and visit 3 as
    the target. Test forecasts use visits 1-2 and the origin-only domain only.
    """
    spacing = tuple(float(s) for s in spacing_mm)
    voxel_mm3 = float(np.prod(spacing))
    # run_outer_patient_cv calls fit once per outer fold in fold order 0..n-1; the count seeds the U-Net per fold.
    fold_counter = {"next": 0}

    def fit(training_records):
        training_records = tuple(training_records)
        fold = fold_counter["next"]
        fold_counter["next"] += 1
        pairs, pair_ids, counts, log_vols, gaps, excluded, outside_crop = [], [], [], [], [], 0, 0
        for r in training_records:
            dom = _mask(domain_by_patient[r.patient_id], name="domain")
            v1, v2, v3 = (_mask(m, name="mask") & dom for m in r.masks[:3])
            gap = float(r.days[2] - r.days[1])
            counts.append(new_label_pieces(v2, v3, spacing_mm=spacing, domain=dom))
            gaps.append(gap)
            log_vols.append(np.log(v2.sum() * voxel_mm3) if v2.any() else np.nan)
            if not v2.any():
                excluded += 1
                continue
            x = unet_inputs(v1, v2, dom)
            if recipe.crop_voxels:
                box = history_crop_box(v2, v2.shape, recipe.crop_voxels)
                outside_crop += int(v3.sum() - v3[box].sum())
                x, v3 = x[(slice(None),) + box], v3[box]
            pairs.append((x, v3))
            pair_ids.append(r.patient_id)
        unet = train_unet_models(pairs, fold=fold, recipe=recipe)
        unet["n_training_empty_history_excluded"] = excluded
        unet["training_fingerprint"] = unet_training_fingerprint(pair_ids, pairs, fold=fold, recipe=recipe)
        unet["training_target_voxels_outside_crop"] = outside_crop
        return {"base": fit_nuisance(training_records), "unet": unet,
                "count": fit_count_model(counts, log_vols, gaps)}

    def _add(prediction, history, fitted):
        dom = _mask(domain_by_patient[history.patient_id], name="domain")
        v1, v2 = (_mask(m, name="history mask") & dom for m in history.masks[:2])
        box = None
        if not v2.any():
            unet_mask, status = v2.copy(), "empty_history_copy_last"
        elif fitted["unet"]["status"] != "ok":
            unet_mask, status = None, "unavailable"
        else:
            x = unet_inputs(v1, v2, dom)
            if recipe.crop_voxels:
                box = history_crop_box(v2, v2.shape, recipe.crop_voxels)
                cropped = predict_unet(fitted["unet"], x[(slice(None),) + box], domain=dom[box])["mask"]
                unet_mask = uncrop_mask(cropped, box, v2.shape)
            else:
                unet_mask = predict_unet(fitted["unet"], x, domain=dom)["mask"]
            status = "forecast"
        log_v = np.log(v2.sum() * voxel_mm3) if v2.any() else np.nan
        if not isinstance(prediction, Mapping):
            prediction = {"fallback_masks": {"hybrid": prediction, "no_pde": prediction}, "history_mask": v2}
        return {**prediction, "_unet_mask": unet_mask, "_unet_status": status, "_unet_crop_box": box,
                "_unet_fingerprint": fitted["unet"].get("training_fingerprint"),
                "_count_rate_per_day": count_rate_per_day(fitted["count"], log_v),
                "_count_mean_rate_per_day": float(fitted["count"]["mean_rate_per_day"]),
                "_count_status": fitted["count"]["status"], "_history_last_day": float(history.days[-1])}

    def wrapped_forecast(history, fitted):
        return _add(forecast(history, fitted["base"]), history, fitted)

    def wrapped_fallback(history, fitted):
        return _add(fallback_forecast(history, fitted["base"]), history, fitted)

    def wrapped_score(prediction, target):
        base = {k: v for k, v in prediction.items() if not k.startswith("_")}
        scores = dict(score(base, target))
        pid = target.patient_id
        dom = _mask(domain_by_patient[pid], name="domain")
        target_mask = _mask(target.mask, name="target")
        support = None if valid_support_by_patient is None else _mask(valid_support_by_patient[pid], name="support")
        scores["unet_status"] = prediction["_unet_status"]
        scores["unet_training_fingerprint"] = prediction["_unet_fingerprint"]
        if prediction["_unet_crop_box"] is not None:
            # Target core outside the history-only crop cannot be predicted; it stays in the target as a miss.
            scores["unet_target_voxels_outside_crop"] = int(
                target_mask.sum() - target_mask[prediction["_unet_crop_box"]].sum())
        if prediction["_unet_mask"] is not None:
            m = prediction["_unet_mask"]
            tm = target_mask
            if support is not None:
                if np.any(m & ~support):
                    raise ValueError("U-Net prediction falls outside validated target support")
                tm = target_mask & support
            scores["dice_unet"] = float(dice_with_empty_rule(m, tm))
            if "dice_hybrid" in scores:
                scores["paired_delta_dice_hybrid_minus_unet"] = float(scores["dice_hybrid"] - scores["dice_unet"])
        # Count model: target day and target mask are used here, after the forecast, never before.
        history_last = _mask(base.get("history_mask"), name="history mask") if base.get("history_mask") is not None \
            else None
        if history_last is None:
            raise ValueError("count scoring needs the visit-2 history mask")
        gap = float(target.day - prediction["_history_last_day"])
        y = new_label_pieces(history_last & dom, target_mask & dom, spacing_mm=spacing, domain=dom)
        mu = prediction["_count_rate_per_day"] * gap
        mu_ref = prediction["_count_mean_rate_per_day"] * gap
        scores.update({"count_observed_new_pieces": int(y), "count_predicted": float(mu),
                       "count_status": prediction["_count_status"], "count_interval_days": gap,
                       **{f"count_{k}": v for k, v in count_metrics(mu, y).items()},
                       **{f"count_mean_rate_{k}": v for k, v in count_metrics(mu_ref, y).items()}})
        return scores

    return fit, wrapped_forecast, wrapped_fallback, wrapped_score
