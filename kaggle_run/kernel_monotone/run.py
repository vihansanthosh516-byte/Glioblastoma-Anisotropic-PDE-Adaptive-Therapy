# Kaggle GPU runner: script 100 production forecasts with the positivity-preserving (Selling) solver (GRAND_PLAN L3).
# Inputs: private dataset vihansanthosh/gbm-forecast-masks-private (built by src/108_pack_masks_for_kaggle.py).
# Stages: forecast (first pairs, full grid) -> select -> selected (later pairs). Analysis runs locally afterwards.
# Each stage runs in shards so a session time-out still leaves finished per-pair caches in /kaggle/working.
import glob
import os
import shutil
import subprocess
import sys
import time

T0 = time.time()
DEADLINE_S = 11.0 * 3600          # Kaggle GPU sessions stop at 12 h
SRC = glob.glob("/kaggle/input/**/mu_masks_1mm.npz", recursive=True)[0]
DS = os.path.dirname(SRC)
REPO = "/kaggle/working/repo"
if not os.path.exists(REPO):
    if os.path.isdir(os.path.join(DS, "repo")):
        shutil.copytree(os.path.join(DS, "repo"), REPO)
    else:   # uploaded with --dir-mode zip and not extracted
        shutil.unpack_archive(os.path.join(DS, "repo.zip"), REPO)
        inner = os.path.join(REPO, "repo")
        if os.path.isdir(inner) and not os.path.exists(os.path.join(REPO, "src")):
            for n in os.listdir(inner):
                shutil.move(os.path.join(inner, n), REPO)
env = dict(os.environ, GBM_MASK_PACK=SRC, GBM_TORCH_DEVICE="cuda", PYTHONUNBUFFERED="1")


def sh(args):
    print(">>", " ".join(args), f"[{(time.time() - T0) / 60:.1f} min]", flush=True)
    r = subprocess.run([sys.executable, "-u", "src/100_pde_manifest.py"] + args, cwd=REPO, env=env)
    if r.returncode != 0:
        raise SystemExit(f"failed: {args}")


subprocess.run(["nvidia-smi"])
import torch  # noqa: E402
print("torch", torch.__version__, "cuda", torch.cuda.is_available(), flush=True)
assert torch.cuda.is_available(), "no GPU: enable the GPU accelerator"

# smoke: one pair with the old solver on GPU (compared locally with the CPU cache)
sh(["--stage", "forecast", "--limit", "1"])

N = 24
for stage in ("forecast", "selected"):
    if stage == "selected":
        sh(["--stage", "select", "--solver", "monotone"])
    for i in range(N):
        if time.time() - T0 > DEADLINE_S:
            print("deadline reached, stopping early", flush=True)
            break
        sh(["--stage", stage, "--solver", "monotone", "--shard", str(i), "--n-shards", str(N)])

shutil.make_archive("/kaggle/working/caches", "zip", os.path.join(REPO, "output"))
print(f"done in {(time.time() - T0) / 3600:.2f} h", flush=True)
