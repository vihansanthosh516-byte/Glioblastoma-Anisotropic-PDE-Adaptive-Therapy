"""Track A region DEVELOPMENT (spent data, exploratory): why does conformal coverage break across sites,
and which recalibration keeps it? Uses GBM-Space donors AT10, AT3, AT4 already spent by the locked test.
Results are development evidence only; any claim needs a fresh locked test.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import trackA_region_transfer as T

OUT = T.PROJECT_ROOT / "output" / "trackA_region_dev_conformal" / "dev_result.json"
ALPHA = 0.10


def lac_scores(p, y):
    return 1 - np.where(y == 1, p, 1 - p)


def qhat(s, w=None, alpha=ALPHA):
    if w is None:
        n = len(s)
        return np.quantile(s, min(1.0, np.ceil((n + 1) * (1 - alpha)) / n), method="higher")
    o = np.argsort(s)
    s, w = s[o], w[o]
    cw = np.cumsum(w) / (w.sum() + w.max())
    k = np.searchsorted(cw, 1 - alpha)
    return s[min(k, len(s) - 1)]


def em_prior(p_t, pi_s, it=200):
    """Saerens et al. 2002: estimate target P(core) from source-calibrated probabilities."""
    pi = pi_s
    for _ in range(it):
        a = p_t * pi / pi_s
        b = (1 - p_t) * (1 - pi) / (1 - pi_s)
        post = a / (a + b)
        new = post.mean()
        if abs(new - pi) < 1e-6:
            break
        pi = new
    return float(pi)


def coverage(p_t, y_t, q1, q0):
    in1 = (1 - p_t) <= q1
    in0 = p_t <= q0
    cov = np.where(y_t == 1, in1, in0).mean()
    return float(cov), float((in1 ^ in0).mean())


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
            s.update(site=site, rep=rep)
            s["y"] = T.label(s.pop("hist"))
            secs[name] = s
        ok = [n for n in secs if (secs[n]["y"] == 1).sum() >= 20 and (secs[n]["y"] == 0).sum() >= 20]
        for src in ok:
            S = secs[src]
            cal = [n for n in ok if n != src and secs[n]["site"] == S["site"]]
            if not cal:
                continue
            ms = S["y"] >= 0
            hv = np.argsort(-S["X"][ms].var(0), kind="stable")[:T.N_HVG]
            hvg = [S["genes"][i] for i in hv]
            Xs, ys = S["X"][ms][:, hv], S["y"][ms]
            model = T.fit_logit(Xs, ys, Xs, rank=True)

            def feats(sec):
                m = sec["y"] >= 0
                gi = {g: i for i, g in enumerate(sec["genes"])}
                A = sec["X"][m]
                Z = np.zeros((A.shape[0], len(hvg)), dtype=A.dtype)
                ok_ = [k for k, g in enumerate(hvg) if g in gi]
                Z[:, ok_] = A[:, [gi[hvg[k]] for k in ok_]]
                return Z, sec["y"][m]

            Xc, yc = feats(secs[cal[0]])
            pc = T.predict(model, Xc)
            sc = lac_scores(pc, yc)
            pi_c = yc.mean()
            for tgt in ok:
                if secs[tgt]["site"] == S["site"]:
                    continue
                Xt, yt = feats(secs[tgt])
                pt = T.predict(model, Xt)
                q = qhat(sc)
                r = {"donor": donor, "source": src, "target": tgt,
                     "true_core_frac_target": float(yt.mean()), "core_frac_cal": float(pi_c)}
                r["marginal"] = coverage(pt, yt, q, q)
                q1, q0 = qhat(sc[yc == 1]), qhat(sc[yc == 0])
                r["class_conditional"] = coverage(pt, yt, q1, q0)
                pi_t = em_prior(pt, pi_c)
                r["em_prior_est"] = pi_t
                w = np.where(yc == 1, pi_t / pi_c, (1 - pi_t) / (1 - pi_c))
                qw = qhat(sc, w)
                r["label_shift_weighted"] = coverage(pt, yt, qw, qw)
                res.append(r)
        del secs
    summ = {}
    for k in ("marginal", "class_conditional", "label_shift_weighted"):
        summ[k] = {"coverage_mean": float(np.mean([r[k][0] for r in res])),
                   "coverage_min": float(np.min([r[k][0] for r in res])),
                   "singleton_mean": float(np.mean([r[k][1] for r in res]))}
    summ["n_pairs"] = len(res)
    summ["em_prior_abs_error_mean"] = float(np.mean([abs(r["em_prior_est"] - r["true_core_frac_target"]) for r in res]))
    summ["core_frac_shift_abs_mean"] = float(np.mean([abs(r["core_frac_cal"] - r["true_core_frac_target"]) for r in res]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"status": "development, spent data", "summary": summ, "pairs": res}, indent=2))
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
