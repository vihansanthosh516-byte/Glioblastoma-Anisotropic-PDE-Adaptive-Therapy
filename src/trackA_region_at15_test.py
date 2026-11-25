"""Track A region: FRESH test of target-adaptive confidence on GBM-Space AT15 (CytAssist v2, untouched).

Protocol: paper/trackA_region_at15_prereg_LOCKED_2026-11-24.md. One run; refuses to overwrite.
Candidate: source-trained rank-normalized classifier (frozen from the locked test), then a
two-component Gaussian mixture on the target section's logits and 6-nearest-neighbour spatial smoothing.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np

import trackA_region_transfer as T
from trackA_region_dev_adaptive import coords, gmm_post, logit, selective, spatial_smooth

OUT = T.PROJECT_ROOT / "output" / "trackA_region_at15_test" / "result.json"
AT15_RE = re.compile(r"^AT15-BRA-4-FO-([A-Z]\d)-S(\d+)\.h5ad$")
FILES = ["AT15-BRA-4-FO-A2-S21.h5ad", "AT15-BRA-4-FO-B2-S25.h5ad", "AT15-BRA-4-FO-C2-S27.h5ad",
         "AT15-BRA-4-FO-C3-S33.h5ad", "AT15-BRA-4-FO-C3-S34.h5ad", "AT15-BRA-4-FO-D3-S36.h5ad",
         "AT15-BRA-4-FO-F1-S18.h5ad", "AT15-BRA-4-FO-F1-S31.h5ad", "AT15-BRA-4-FO-G1-S21.h5ad"]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def has_labels(path):
    f = h5py.File(path, "r")
    ft = T._strs(f["var"]["feature_types"])
    hi = np.where(ft == "Histopath annotation overlap")[0]
    X = f["X"]
    if not len(hi):
        return False
    if isinstance(X, h5py.Group):
        return bool(np.isin(X["indices"][:], hi).any())
    return bool(np.nansum(X[:, hi.min():hi.max() + 1]) > 0)


def matched(score, y, frac):
    k = max(1, int(round(frac * len(y))))
    top = np.argsort(-np.abs(score - 0.5))[:k]
    return float(((score[top] >= 0.5).astype(int) == y[top]).mean()), k / len(y)


def main():
    if OUT.exists():
        print("refusing to overwrite", OUT)
        return 1
    secs = {}
    for name in FILES:
        p = T.DATA / name
        if not has_labels(p):
            continue
        m = AT15_RE.match(name)
        s = T.load_section(p)
        s.update(site=m.group(1), xy=coords(p, s["X"].shape[0]))
        s["y"] = T.label(s.pop("hist"))
        secs[name] = s
    ok = [n for n in secs if (secs[n]["y"] == 1).sum() >= T.MIN_CLASS and (secs[n]["y"] == 0).sum() >= T.MIN_CLASS]
    res = []
    for src in ok:
        S = secs[src]
        ms = S["y"] >= 0
        hv = np.argsort(-S["X"][ms].var(0), kind="stable")[:T.N_HVG]
        hvg = [S["genes"][i] for i in hv]
        Xs = S["X"][ms][:, hv]
        model = T.fit_logit(Xs, S["y"][ms], Xs, rank=True)
        for tgt in ok:
            Tt = secs[tgt]
            if Tt["site"] == S["site"]:
                continue
            gi = {g: i for i, g in enumerate(Tt["genes"])}
            Z = np.zeros((Tt["X"].shape[0], len(hvg)), dtype=Tt["X"].dtype)
            k_ok = [k for k, g in enumerate(hvg) if g in gi]
            Z[:, k_ok] = Tt["X"][:, [gi[hvg[k]] for k in k_ok]]
            p_all = T.predict(model, Z)
            m = Tt["y"] >= 0
            y = Tt["y"][m]
            g_all = gmm_post(logit(p_all))
            cand = selective(spatial_smooth(g_all, Tt["xy"])[m], y)
            r = {"source": src, "target": tgt, "n_scored": int(m.sum()),
                 "candidate": cand,
                 "abl_no_spatial": selective(g_all[m], y),
                 "abl_no_gmm": selective(spatial_smooth(p_all, Tt["xy"])[m], y),
                 "source_prob": selective(p_all[m], y),
                 "matched_source_prob": matched(p_all[m], y, cand[1]),
                 "matched_source_prob_spatial": matched(spatial_smooth(p_all, Tt["xy"])[m], y, cand[1])}
            res.append(r)
    out = {"labelled_sections": sorted(secs), "eligible_sections": ok, "n_pairs": len(res), "pairs": res}
    if res:
        mean = lambda k: float(np.nanmean([r[k][0] for r in res]))
        cm = mean("candidate")
        b1, b2 = mean("matched_source_prob"), mean("matched_source_prob_spatial")
        w1 = np.mean([r["candidate"][0] > r["matched_source_prob"][0] for r in res])
        w2 = np.mean([r["candidate"][0] > r["matched_source_prob_spatial"][0] for r in res])
        out.update({
            "candidate_acc_confident_mean": cm,
            "candidate_frac_confident_mean": float(np.mean([r["candidate"][1] for r in res])),
            "matched_source_prob_mean": b1, "matched_source_prob_spatial_mean": b2,
            "abl_no_spatial_mean": mean("abl_no_spatial"), "abl_no_gmm_mean": mean("abl_no_gmm"),
            "source_prob_mean": mean("source_prob"),
            "win_rate_vs_matched_source_prob": float(w1), "win_rate_vs_matched_source_prob_spatial": float(w2),
        })
        out["verdict"] = "PASS" if (len(res) >= 6 and cm >= 0.90 and cm > b1 and cm > b2
                                    and w1 >= 2 / 3 and w2 >= 2 / 3) else "FAIL"
    else:
        out["verdict"] = "INCONCLUSIVE (no eligible cross-block pairs)"
    out["input_sha256"] = {n: sha(T.DATA / n) for n in secs}
    out["code_sha256"] = {Path(f).name: sha(Path(f)) for f in [__file__, T.__file__,
                          str(T.PROJECT_ROOT / "src" / "trackA_region_dev_adaptive.py")]}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps({k: v for k, v in out.items() if k not in ("pairs", "input_sha256")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
