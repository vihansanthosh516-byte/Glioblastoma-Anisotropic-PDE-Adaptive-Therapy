# Fetch the published PREDICT-GBM U-Net weights (MIT, huggingface.co/LZimmer/PREDICT-GBM-Models) inside Kaggle only.
# Untrusted file: this job only downloads, extracts and LISTS contents. Nothing is unpickled or executed here.
import os
import subprocess
import tarfile
import zipfile

URL = "https://huggingface.co/LZimmer/PREDICT-GBM-Models/resolve/main/unet_weights.tar.gz"
os.makedirs("/kaggle/working/unet_weights", exist_ok=True)
subprocess.run(["wget", "-q", "-O", "/kaggle/tmp_w.tar.gz", URL], check=True)
print("size MB", os.path.getsize("/kaggle/tmp_w.tar.gz") / 1e6)
with tarfile.open("/kaggle/tmp_w.tar.gz") as t:
    members = t.getmembers()
    for m in members:
        print(("DIR " if m.isdir() else "FILE"), m.size, m.name)
    safe = [m for m in members if (m.isfile() or m.isdir()) and not m.name.startswith("/") and ".." not in m.name]
    t.extractall("/kaggle/working/unet_weights", members=safe, filter="data")
for root, _, files in os.walk("/kaggle/working/unet_weights"):
    for f in files:
        p = os.path.join(root, f)
        if zipfile.is_zipfile(p):   # torch checkpoints are zip archives: list entries without unpickling
            with zipfile.ZipFile(p) as z:
                names = z.namelist()
                print("ZIP", p, len(names), names[:12])
        elif f.endswith((".py", ".json", ".yaml", ".yml", ".txt", ".md")):
            print("TEXT", p)
            print(open(p, errors="replace").read()[:3000])
os.remove("/kaggle/tmp_w.tar.gz")
