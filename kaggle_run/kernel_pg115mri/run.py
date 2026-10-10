# Kaggle GPU runner: script 115 (H-2 learned U-Net and hybrid arms, PREDICT-GBM development, 5-fold out-of-fold).
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
code = glob.glob("/kaggle/input/**/src/115_predictgbm_unet_dev.py", recursive=True)
if code:
    shutil.copytree(os.path.dirname(os.path.dirname(code[0])), REPO, dirs_exist_ok=True)
else:
    copy_repo(os.path.dirname(glob.glob("/kaggle/input/**/repo*.zip", recursive=True)[-1]))
assert os.path.exists(os.path.join(REPO, "src", "115_predictgbm_unet_dev.py"))
env = dict(os.environ, GBM_PG_PACK=PACK, GBM_PG_TEST_IDS=os.path.join(DS, "test_ids.json"), GBM_PG_MRI=glob.glob("/kaggle/input/**/pg_mri2mm.npz", recursive=True)[0], PYTHONUNBUFFERED="1")
subprocess.run(["nvidia-smi"])
r = subprocess.run([sys.executable, "-u", "src/115_predictgbm_unet_dev.py", "--arms", "unet_mri,hybrid_mri", "--epochs", "60", "--merge-prev"],
                   cwd=REPO, env=env)
for f in ("unet_dev_mri.json", "unet_dev_mri_pairs.csv"):
    p = os.path.join(REPO, "output", "predictgbm", f)
    if os.path.exists(p):
        shutil.copy2(p, "/kaggle/working/" + f)
print(f"exit {r.returncode}; done in {(time.time() - T0) / 3600:.2f} h", flush=True)
