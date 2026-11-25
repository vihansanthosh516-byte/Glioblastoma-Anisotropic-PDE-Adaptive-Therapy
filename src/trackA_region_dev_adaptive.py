"""Track A region DEVELOPMENT (spent data): target-adaptive confidence for one-section transfer.

Goal: among spots the method calls "confident", accuracy should reach >= 90% on a new site,
while keeping as many spots as possible. Spent donors AT10, AT3, AT4 only. Exploratory.
"""
from __future__ import annotations

import json
import sys

import h5py
import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors

import trackA_region_transfer as T

OUT = T.PROJECT_ROOT / "output" / "trackA_region_dev_adaptive" / "dev_result.json"
CONF = 0.90
K_SP = 6


def coords(path, n_expected):
    f = h5py.File(path, "r")
    c = f["obsm"]["spatial"][:]
    assert c.shape[0] == n_expected
    return c.astype(float)


def logit(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def gmm_post(z):
    """Two-component 1-D Gaussian mixture on target logits; returns P(core) using the higher-mean component."""
    g = GaussianMixture(2, random_state=0, n_init=3).fit(z[:, None])
    hi = int(np.argmax(g.means_.ravel()))
    return g.predict_proba(z[:, None])[:, hi]


def spatial_smooth(p, xy):
    nn = NearestNeighbors(n_neighbors=K_SP + 1).fit(xy)
    _, idx = nn.kneighbors(xy)
    return p[idx].mean(1)


def selective(post, y, conf=CONF):
    """Confident if max(post, 1-post) >= conf. Returns (accuracy among confident, fraction confident)."""
    c = np.maximum(post, 1 - post) >= conf
    if c.sum() == 0:
        return float("nan"), 0.0
    pred = (post >= 0.5).astype(int)
    return float((pred[c] == y[c]).mean()), float(c.mean())


def conformal_single(p_cal, y_cal, p_t, y_t, alpha=0.10):
    s = 1 - np.where(y_cal == 1, p_cal, 1 - p_cal)
    n = len(s)
    q = np.quantile(s, min(1.0, np.ceil((n + 1) * (1 - alpha)) / n), method="higher")
    in1, in0 = (1 - p_t) <= q, p_t <= q
    single = in1 ^ in0
    if single.sum() == 0:
        return float("nan"), 0.0
    pred = in1.astype(int)
    return float((pred[single] == y_t[single]).mean()), float(single.mean())


def main():
    rows = [l.rstrip("\n").split("\t") for l in open(T.COVER)]
    by = {}
    for _d, name, has, _n in rows:
        m = T.NAME_RE.match(name)
        if has == "True" and m and m.group(1) in ("AT10", "AT3", "AT4"):
            by.setdefault(m.group(1), []).append((name, int(m.group(2)), int(m.group(3))))
    res = []
    for donor, items in sorted(by.items()):
        secs = {}
        for name, site, rep in items:
            s = T.load_section(T.DATA / name)
            s.update(site=site, rep=rep, xy=coords(T.DATA / name, s["X"].shape[0]))
            s["y"] = T.label(s.pop("hist"))
            secs[name] = s
        ok = [n for n in secs if (secs[n]["y"] == 1).sum() >= 20 and (secs[n]["y"] == 0).sum() >= 20]
        for src in ok:
            S = secs[src]
            cal = [n for n in ok if n != src and secs[n]["site"] == S["site"]]
            ms = S["y"] >= 0
            hv = np.argsort(-S["X"][ms].var(0), kind="stable")[:T.N_HVG]
            hvg = [S["genes"][i] for i in hv]
            Xs_ = S["X"][ms][:, hv]
            model = T.fit_logit(Xs_, S["y"][ms], Xs_, rank=True)

            def feats(sec, all_spots=False):
                gi = {g: i for i, g in enumerate(sec["genes"])}
                A = sec["X"] if all_spots else sec["X"][sec["y"] >= 0]
                Z = np.zeros((A.shape[0], len(hvg)), dtype=A.dtype)
                k_ok = [k for k, g in enumerate(hvg) if g in gi]
                Z[:, k_ok] = A[:, [gi[hvg[k]] for k in k_ok]]
                return Z

            for tgt in ok:
                Tt = secs[tgt]
                if Tt["site"] == S["site"]:
                    continue
                # predict on ALL target spots (unlabeled use is allowed), score on labelled ones
                p_all = T.predict(model, feats(Tt, all_spots=True))
                m = Tt["y"] >= 0
                y = Tt["y"][m]
                r = {"donor": donor, "source": src, "target": tgt}
                r["source_prob"] = selective(p_all[m], y)
                g_all = gmm_post(logit(p_all))
                r["target_gmm"] = selective(g_all[m], y)
                r["target_gmm_spatial"] = selective(spatial_smooth(g_all, Tt["xy"])[m], y)
                r["source_prob_spatial"] = selective(spatial_smooth(p_all, Tt["xy"])[m], y)
                # matched-coverage baselines: same fraction of confident spots as the candidate
                gs = spatial_smooth(g_all, Tt["xy"])[m]
                frac = r["target_gmm_spatial"][1]
                for nm, sc in (("source_prob_matched", p_all[m]), ("source_prob_spatial_matched", spatial_smooth(p_all, Tt["xy"])[m])):
                    k = max(1, int(round(frac * len(y))))
                    top = np.argsort(-np.abs(sc - 0.5))[:k]
                    r[nm] = (float(((sc[top] >= 0.5).astype(int) == y[top]).mean()), k / len(y))
                if cal:
                    C = secs[cal[0]]
                    r["replicate_conformal"] = conformal_single(T.predict(model, feats(C)), C["y"][C["y"] >= 0], p_all[m], y)
                res.append(r)
        del secs
    arms = ["source_prob", "source_prob_spatial", "replicate_conformal", "target_gmm", "target_gmm_spatial", "source_prob_matched", "source_prob_spatial_matched"]
    summ = {}
    for a in arms:
        v = [r[a] for r in res if a in r and not np.isnan(r[a][0])]
        summ[a] = {"acc_confident_mean": float(np.mean([x[0] for x in v])),
                   "acc_confident_min": float(np.min([x[0] for x in v])),
                   "pairs_acc_ge_0.90": int(sum(x[0] >= 0.90 for x in v)),
                   "frac_confident_mean": float(np.mean([x[1] for x in v])),
                   "n_pairs": len(v)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"status": "development, spent data", "summary": summ, "pairs": res}, indent=2))
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    sys.exit(main())
