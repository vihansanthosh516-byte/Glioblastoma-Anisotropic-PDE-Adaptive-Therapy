"""Rebuild of the released PREDICT-GBM U-Net (MONAI UNet layout, num_res_units=0) without MONAI.

State-dict keys match MONAI: model.0 (down), model.1 (SkipConnection.submodule = Sequential(down, skip..., up)), model.2
(final transposed conv, conv only). Channels 32-64-128-256-512-1024, strides 2, kernel 3. The checkpoint has no norm or
activation parameters, so norm is InstanceNorm3d(affine=False) or none and the activation is parameter-free; the variant
is chosen by matching the released out-of-fold predictions on development patients (script 125 --stage validate).
Weights are loaded with torch.load(weights_only=True) on Kaggle only.
"""
from __future__ import annotations

import torch
from torch import nn

CH = (32, 64, 128, 256, 512, 1024)


class Conv(nn.Module):
    def __init__(self, cin, cout, stride, transposed, conv_only, norm, act):
        super().__init__()
        if transposed:
            self.conv = nn.ConvTranspose3d(cin, cout, 3, stride, padding=1, output_padding=stride - 1)
        else:
            self.conv = nn.Conv3d(cin, cout, 3, stride, padding=1)
        layers = []
        if not conv_only:
            if norm == "instance":
                layers.append(nn.InstanceNorm3d(cout, affine=False))
            layers.append({"relu": nn.ReLU(), "leaky": nn.LeakyReLU(0.01), "prelu0": nn.LeakyReLU(0.25)}[act])
        self.adn = nn.Sequential(*layers)

    def forward(self, x):
        return self.adn(self.conv(x))


class Skip(nn.Module):
    def __init__(self, sub):
        super().__init__()
        self.submodule = sub

    def forward(self, x):
        return torch.cat([x, self.submodule(x)], 1)


def _level(cin, chs, norm, act, top):
    c = chs[0]
    if len(chs) == 2:   # bottom: down conv then bottom conv (stride 1)
        sub = Conv(c, chs[1], 1, False, False, norm, act)
        up_in = c + chs[1]
    else:
        sub = _level(c, chs[1:], norm, act, False)
        up_in = c * 2
    down = Conv(cin, c, 2, False, False, norm, act)
    up = Conv(up_in, OUT if top else cin, 2, True, top, norm, act)
    return nn.Sequential(down, Skip(sub), up)


OUT = 1


class PublishedUNet(nn.Module):
    def __init__(self, cin=7, norm="instance", act="relu"):
        super().__init__()
        self.model = _level(cin, CH, norm, act, True)

    def forward(self, x):
        return self.model(x)
