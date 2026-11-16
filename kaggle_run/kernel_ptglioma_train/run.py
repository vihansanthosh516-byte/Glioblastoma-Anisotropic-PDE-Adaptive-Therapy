# Track B lock chain step 3: train ONE two-sequence (T1Gd, FLAIR) nnU-Net on MU + UCSF masks (spent cohorts only).
# Protocol: paper/trackB_vig_cfb_waw_protocol_2026-11-15.md section 4. 3d_fullres, fold "all", 250 epochs, last checkpoint, no model selection.
# One session cannot hold 250 epochs. The run stops itself before the session limit and a later version resumes from checkpoint_latest.pth
# (add this kernel to kernel_sources of the next version). No test-cohort file is read here.
import glob, hashlib, json, os, shutil, subprocess, sys, time, zipfile
T0 = time.time(); BUDGET_H = 10.3
W = "/kaggle/working"; PRE = "/tmp/nnunet_pre"; RES = f"{W}/nnUNet_results"
run = lambda c, **k: subprocess.run(c, shell=isinstance(c, str), **k)
run("nvidia-smi")
# nnunetv2 2.5.1 does not run with the torch of the Kaggle image (2.11): its LR scheduler call fails. Use the current release; the
# installed version is written to train_status.json and must be the same in every later session and at inference.
NNUNET_PIN = os.environ.get("NNUNET_VERSION", "2.8.1")           # version of session 1 (train_status.json of kernel version 3)
run(f"{sys.executable} -m pip install -q nnunetv2=={NNUNET_PIN}", check=True)
import nnunetv2, torch
ver = {"nnunetv2": getattr(nnunetv2, "__version__", None) or subprocess.run(f"{sys.executable} -m pip show nnunetv2", shell=True, capture_output=True, text=True).stdout.split("Version:")[1].split()[0],
       "torch": torch.__version__, "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}
print("versions", ver, flush=True)

hit = glob.glob("/kaggle/input/**/Dataset501_PTGlioma/dataset.json", recursive=True)
if not hit:                                                   # directory uploaded as a zip that was not unpacked
    for z in glob.glob("/kaggle/input/**/*.zip", recursive=True):
        zipfile.ZipFile(z).extractall("/tmp/raw_unzip")
    hit = glob.glob("/tmp/raw_unzip/**/Dataset501_PTGlioma/dataset.json", recursive=True)
SRC = os.path.dirname(hit[0])
# Kaggle may unpack .nii.gz to .nii. nnU-Net needs one file ending, so build a raw folder of links with a matching dataset.json.
end = ".nii.gz" if glob.glob(f"{SRC}/labelsTr/*.nii.gz") else ".nii"
RAW = "/tmp/nnunet_raw"; D = f"{RAW}/Dataset501_PTGlioma"; os.makedirs(D, exist_ok=True)
for sub in ("imagesTr", "labelsTr"):
    if not os.path.exists(f"{D}/{sub}"): os.symlink(f"{SRC}/{sub}", f"{D}/{sub}")
dj = json.load(open(hit[0])); dj["file_ending"] = end; json.dump(dj, open(f"{D}/dataset.json", "w"))
n_img, n_lab = len(glob.glob(f"{D}/imagesTr/*{end}")), len(glob.glob(f"{D}/labelsTr/*{end}"))
print("raw", SRC, "ending", end, "images", n_img, "labels", n_lab, flush=True); assert n_img == 2 * n_lab and n_lab == dj["numTraining"]

old = glob.glob("/kaggle/input/**/nnUNet_results", recursive=True)   # checkpoint of an earlier session, if this version has one as input
if old and not os.path.exists(RES):
    shutil.copytree(old[0], RES); run(f"chmod -R u+w {RES}"); print("restored results from", old[0], flush=True)
os.makedirs(PRE, exist_ok=True); os.makedirs(RES, exist_ok=True)
env = dict(os.environ, nnUNet_raw=RAW, nnUNet_preprocessed=PRE, nnUNet_results=RES, nnUNet_n_proc_DA="4", nnUNet_compile="f")

tr_dir = os.path.join(os.path.dirname(nnunetv2.__file__), "training", "nnUNetTrainer", "variants")
open(os.path.join(tr_dir, "nnUNetTrainer_VIG250.py"), "w").write('''import random, numpy as np, torch
from nnunetv2.training.nnUNetTrainer.nnUNetTrainer import nnUNetTrainer
random.seed(20261115); np.random.seed(20261115); torch.manual_seed(20261115)


class nnUNetTrainer_VIG250(nnUNetTrainer):
    """Default trainer. Only changes: 250 epochs, checkpoint every 5 epochs, no end-of-training validation pass (fold all has no held-out cases)."""
    def initialize(self):
        self.num_epochs = 250
        self.save_every = 5
        super().initialize()

    def perform_actual_validation(self, *a, **k):
        return None
''')
run("nnUNetv2_plan_and_preprocess -d 501 -c 3d_fullres -np 4 -npfp 4", env=env, check=True)
print("preprocess done after h", round((time.time() - T0) / 3600, 2), flush=True)
shutil.copy(glob.glob(f"{PRE}/Dataset501_PTGlioma/nnUNetPlans.json")[0], f"{W}/nnUNetPlans.json")

left = BUDGET_H * 3600 - (time.time() - T0); status = "finished"
log = open(f"{W}/train_stdout.txt", "a")
p = subprocess.Popen("nnUNetv2_train 501 3d_fullres all -tr nnUNetTrainer_VIG250 --c", shell=True, env=env, stdout=log, stderr=subprocess.STDOUT)
try:
    rc = p.wait(timeout=max(left, 60)); status = "finished" if rc == 0 else f"error_rc_{rc}"
except subprocess.TimeoutExpired:
    p.kill(); run("pkill -f nnUNetv2_train"); status = "stopped_for_session_limit"
log.close(); run(f"tail -n 30 {W}/train_stdout.txt")

fold = glob.glob(f"{RES}/Dataset501_PTGlioma/*/fold_all"); out = {"status": status, "versions": ver, "raw_cases": n_lab, "hours": round((time.time() - T0) / 3600, 2)}
if fold:
    logs = sorted(glob.glob(f"{fold[0]}/training_log_*.txt")); ep = [l for f in logs for l in open(f) if "Epoch " in l and "time" not in l]
    out["last_epoch_line"] = ep[-1].strip() if ep else None
    for c in ("checkpoint_final.pth", "checkpoint_latest.pth"):
        f = f"{fold[0]}/{c}"
        if os.path.exists(f): out[c + "_sha256"] = hashlib.sha256(open(f, "rb").read()).hexdigest()
    for junk in glob.glob(f"{fold[0]}/checkpoint_best.pth"): os.remove(junk)      # protocol uses the last checkpoint only
json.dump(out, open(f"{W}/train_status.json", "w"), indent=1); print(json.dumps(out, indent=1))
