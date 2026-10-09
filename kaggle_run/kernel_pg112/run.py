# Kaggle GPU runner: script 112 (H-2 models M1/M2 on PREDICT-GBM development patients, exact solver).
# Input: private dataset vihansanthosh/predictgbm-inputs-private (src/113_pack_predictgbm_for_kaggle.py).
# Test recurrences are not in the pack; the loader lock refuses them anyway.
import glob
import os
import shutil
import subprocess
import sys
import time

T0 = time.time()
DEADLINE_S = 11.0 * 3600
PACK = glob.glob("/kaggle/input/**/pg_inputs.npz", recursive=True)[0]
DS = os.path.dirname(PACK)
REPO = "/kaggle/working/repo"
if not os.path.exists(REPO):
    if os.path.isdir(os.path.join(DS, "repo")):
        shutil.copytree(os.path.join(DS, "repo"), REPO)
    else:
        shutil.unpack_archive(os.path.join(DS, "repo.zip"), REPO)
        inner = os.path.join(REPO, "repo")
        if os.path.isdir(inner) and not os.path.exists(os.path.join(REPO, "src")):
            for n in os.listdir(inner):
                shutil.move(os.path.join(inner, n), REPO)
env = dict(os.environ, GBM_PG_PACK=PACK, GBM_PG_TEST_IDS=os.path.join(DS, "test_ids.json"),
           GBM_TORCH_DEVICE="cuda", PYTHONUNBUFFERED="1")
subprocess.run(["nvidia-smi"])
N = 24
for i in range(N):
    if time.time() - T0 > DEADLINE_S:
        print("deadline reached", flush=True)
        break
    print(f">> shard {i} [{(time.time() - T0) / 60:.1f} min]", flush=True)
    r = subprocess.run([sys.executable, "-u", "src/112_predictgbm_models_dev.py", "--stage", "run", "--shard", str(i),
                        "--n-shards", str(N)], cwd=REPO, env=env)
    if r.returncode != 0:
        raise SystemExit(f"shard {i} failed")
shutil.make_archive("/kaggle/working/cache_112", "zip", os.path.join(REPO, "output", "predictgbm"))
print(f"done in {(time.time() - T0) / 3600:.2f} h", flush=True)
