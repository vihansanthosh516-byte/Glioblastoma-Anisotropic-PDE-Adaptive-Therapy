# Kaggle GPU runner: script 114 (H-2 M3 per-patient fit, PREDICT-GBM development).
# Inputs: private datasets predictgbm-inputs-private (data, atlas; no test recurrences) and gbm-code-private (code only).
import glob
import os
import shutil
import subprocess
import sys
import time

T0 = time.time()
PACK = glob.glob("/kaggle/input/**/pg_inputs.npz", recursive=True)[0]
DS = os.path.dirname(PACK)
REPO = "/kaggle/working/repo"


def copy_repo(src_dir):
    if os.path.isdir(os.path.join(src_dir, "repo")):
        shutil.copytree(os.path.join(src_dir, "repo"), REPO, dirs_exist_ok=True)
    else:
        tmp = "/kaggle/working/_tmp"
        shutil.unpack_archive(glob.glob(os.path.join(src_dir, "repo*.zip"))[0], tmp)
        inner = os.path.join(tmp, "repo") if os.path.isdir(os.path.join(tmp, "repo")) else tmp
        shutil.copytree(inner, REPO, dirs_exist_ok=True)
        shutil.rmtree(tmp)


copy_repo(DS)
code = glob.glob("/kaggle/input/**/src/114_predictgbm_m3_patient_fit.py", recursive=True)
if code:
    shutil.copytree(os.path.dirname(os.path.dirname(code[0])), REPO, dirs_exist_ok=True)
else:
    copy_repo(os.path.dirname(glob.glob("/kaggle/input/**/repo*.zip", recursive=True)[-1]))
assert os.path.exists(os.path.join(REPO, "src", "114_predictgbm_m3_patient_fit.py"))
env = dict(os.environ, GBM_PG_PACK=PACK, GBM_PG_TEST_IDS=os.path.join(DS, "test_ids.json"), GBM_TORCH_DEVICE="cuda", PYTHONUNBUFFERED="1")
subprocess.run(["nvidia-smi"])
r = None
for i in range(24):
    print(f">> shard {i} [{(time.time() - T0) / 60:.1f} min]", flush=True)
    r = subprocess.run([sys.executable, "-u", "src/114_predictgbm_m3_patient_fit.py", "--stage", "run", "--shard", str(i), "--n-shards", "24"], cwd=REPO, env=env)
    if r.returncode != 0:
        raise SystemExit(f"shard {i} failed")
shutil.make_archive("/kaggle/working/cache_114", "zip", os.path.join(REPO, "output", "predictgbm"), "cache_114")
print(f"exit {r.returncode}; done in {(time.time() - T0) / 3600:.2f} h", flush=True)
