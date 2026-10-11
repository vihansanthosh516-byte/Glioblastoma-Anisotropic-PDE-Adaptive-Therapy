# Inspect the published U-Net checkpoints with torch.load(weights_only=True): tensors and plain containers only,
# no arbitrary code. Prints top-level keys, parameter names and shapes of one checkpoint.
import glob
import torch

fs = sorted(glob.glob("/kaggle/input/**/ens_1_crv_1.tar", recursive=True))
print(fs)
ck = torch.load(fs[0], map_location="cpu", weights_only=True)
print(type(ck))
if isinstance(ck, dict):
    for k, v in ck.items():
        print("TOP", k, type(v).__name__, (len(v) if hasattr(v, "__len__") else v) if not torch.is_tensor(v) else tuple(v.shape))
    sd = ck.get("model_state_dict") or ck.get("state_dict") or ck.get("model") or ck
    n = 0
    for k, v in sd.items():
        if torch.is_tensor(v):
            n += v.numel()
            print(k, tuple(v.shape))
    print("params", n)
    for k, v in ck.items():
        if not torch.is_tensor(v) and not isinstance(v, dict):
            print("META", k, str(v)[:2000])
