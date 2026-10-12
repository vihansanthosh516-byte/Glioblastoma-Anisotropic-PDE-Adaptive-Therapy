# Script 125 stage validate (CPU): pick the published U-Net rebuild variant that reproduces the RELEASED out-of-fold
# predictions on 8 TUM dev patients (no recurrence used). Weights from kernel gbm-unet-weights-fetch, weights_only=True.
import glob, itertools, json, sys, time
import numpy as np
import torch
from scipy.stats import spearmanr
sys.path.insert(0, glob.glob("/kaggle/input/**/unet_published.py", recursive=True)[0].rsplit("/", 1)[0])
from unet_published import PublishedUNet
torch.set_num_threads(4)
V = np.load(glob.glob("/kaggle/input/**/unetval.npz", recursive=True)[0])
P = np.load(glob.glob("/kaggle/input/**/pg_inputs.npz", recursive=True)[0])
W = {f.rsplit("/", 1)[1][:-4]: f for f in glob.glob("/kaggle/input/**/ens_*_crv_*.tar", recursive=True)}
print(sorted(W))
ids = sorted({k.split("|")[0] for k in V.files})


def bits(pid, n):
    k = f"{pid}|{n}"
    if k + "|bits" in P.files:
        return np.unpackbits(P[k + "|bits"])[: 240 * 240 * 155].reshape(240, 240, 155)
    return P[k]


def inputs(pid, segmode):
    seg = np.rint(bits(pid, "tumor_seg.nii.gz")).astype(np.float32)
    s = {"raw": seg, "div3": seg / 3.0, "div4": seg / 4.0}[segmode]
    ch = [s] + [bits(pid, f"{t}_pbmap.nii.gz").astype(np.float32) for t in ("wm", "gm", "csf")] + \
         [V[f"{pid}|t1c_bet_normalized"].astype(np.float32), V[f"{pid}|flair_bet_normalized"].astype(np.float32),
          np.zeros((240, 240, 155), np.float32)]
    x = np.zeros((7, 256, 256, 160), np.float32)
    x[:, 8:248, 8:248, 2:157] = np.stack(ch)
    return torch.from_numpy(x)[None]


def predict(net, x):
    with torch.no_grad():
        return net(x)[0, 0, 8:248, 8:248, 2:157].numpy()


def load(name, norm, act):
    net = PublishedUNet(7, norm, act)
    net.load_state_dict(torch.load(W[name], map_location="cpu", weights_only=True))
    return net.eval()


def score(pred, pid):
    ref = V[f"{pid}|unet_pred"].astype(np.float32)
    m = bits(pid, "t1c_bet_mask.nii.gz") > 0
    return float(spearmanr(pred[m][::50], ref[m][::50])[0]), float(np.abs(pred[m] - ref[m]).mean())


t0 = time.time()
res = {}
for norm, act, segmode in itertools.product(("instance", "none"), ("relu", "leaky", "prelu0"), ("raw", "div3", "div4")):
    r = []
    for pid in ids[:2]:
        k = int(V[f"{pid}|fold"]) + 1
        r.append(score(predict(load(f"ens_1_crv_{k}", norm, act), inputs(pid, segmode)), pid))
    res[f"{norm}|{act}|{segmode}"] = r
    print(norm, act, segmode, [(round(a, 3), round(b, 3)) for a, b in r], f"{time.time() - t0:.0f}s", flush=True)
best = max(res, key=lambda v: np.mean([a for a, _ in res[v]]))
norm, act, segmode = best.split("|")
print("BEST", best, flush=True)
full = {}
for pid in ids:   # confirm: mean of the 3 members of the patient's held-out fold, all 8 patients
    k = int(V[f"{pid}|fold"]) + 1
    x = inputs(pid, segmode)
    pred = np.mean([predict(load(f"ens_{e}_crv_{k}", norm, act), x) for e in (1, 2, 3)], 0)
    full[pid] = score(pred, pid)
    print(pid, k, full[pid], flush=True)
json.dump({"screen": res, "best": best, "confirm": full}, open("/kaggle/working/unetval.json", "w"), indent=1)
print(f"done {(time.time() - t0) / 60:.1f} min")
