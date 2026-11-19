# Track B lock chain step 3b: run the frozen two-sequence nnU-Net on UCSD scans (spent cohort, not in training) for gate G1.
# Weights: checkpoint_final.pth of kernel trackb-ptglioma-nnunet-train version 4 (nnunetv2 2.8.1). No test-cohort file is read.
import glob, hashlib, json, os, shutil, subprocess, sys, zipfile
run = lambda c, **k: subprocess.run(c, shell=isinstance(c, str), **k)
run("nvidia-smi")
run(f"{sys.executable} -m pip install -q nnunetv2==2.8.1", check=True)
import nnunetv2, torch
W = "/kaggle/working"; RES = "/tmp/nnunet_results"; OUT = f"{W}/pred"; os.makedirs(OUT, exist_ok=True)
for z in glob.glob("/kaggle/input/**/*.zip", recursive=True):     # in case a directory was kept as a zip
    zipfile.ZipFile(z).extractall("/tmp/unz")
roots = ["/kaggle/input", "/tmp/unz"]
def find(pat):
    for r in roots:
        h = glob.glob(f"{r}/**/{pat}", recursive=True)
        if h: return h
    return []
ck = find("checkpoint_final.pth")[0]; md = os.path.dirname(ck)
imgs = find("*_0000.nii*"); IMD = os.path.dirname(imgs[0]); end = ".nii.gz" if imgs[0].endswith(".gz") else ".nii"
n0, n1 = len(glob.glob(f"{IMD}/*_0000{end}")), len(glob.glob(f"{IMD}/*_0001{end}")); print("images", IMD, end, n0, n1, flush=True); assert n0 == n1 > 0
sha = hashlib.sha256(open(ck, "rb").read()).hexdigest(); print("checkpoint sha256", sha, flush=True)
M = f"{RES}/Dataset501_PTGlioma/nnUNetTrainer_VIG250__nnUNetPlans__3d_fullres"; os.makedirs(f"{M}/fold_all", exist_ok=True)
for f in ("dataset.json", "dataset_fingerprint.json", "plans.json"): shutil.copyfile(f"{md}/{f}", f"{M}/{f}")
os.symlink(ck, f"{M}/fold_all/checkpoint_final.pth")
dj = json.load(open(f"{M}/dataset.json")); dj["file_ending"] = end; json.dump(dj, open(f"{M}/dataset.json", "w"))
tr = os.path.join(os.path.dirname(nnunetv2.__file__), "training", "nnUNetTrainer", "variants")
open(os.path.join(tr, "nnUNetTrainer_VIG250.py"), "w").write('''from nnunetv2.training.nnUNetTrainer.nnUNetTrainer import nnUNetTrainer


class nnUNetTrainer_VIG250(nnUNetTrainer):
    def initialize(self):
        self.num_epochs = 250
        self.save_every = 5
        super().initialize()
''')
env = dict(os.environ, nnUNet_results=RES, nnUNet_raw="/tmp/x", nnUNet_preprocessed="/tmp/y")
rc = run(f"nnUNetv2_predict -i {IMD} -o {OUT} -d 501 -c 3d_fullres -tr nnUNetTrainer_VIG250 -p nnUNetPlans -f all -chk checkpoint_final.pth -device cuda", env=env).returncode
preds = glob.glob(f"{OUT}/*{end}")
json.dump({"returncode": rc, "n_inputs": n0, "n_predictions": len(preds), "checkpoint_sha256": sha, "nnunetv2": nnunetv2.__version__ if hasattr(nnunetv2, "__version__") else "2.8.1",
           "torch": torch.__version__}, open(f"{W}/g1_status.json", "w"), indent=1)
print(open(f"{W}/g1_status.json").read())
