import torch
import torch.nn as nn
import torch.nn.functional as F


class SpectralConv3d(nn.Module):
    """3D Fourier Layer: FFT -> Linear Transform on lower modes -> Inverse FFT"""
    def __init__(self, in_channels, out_channels, modes1, modes2, modes3):
        super(SpectralConv3d, self).__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.modes1 = modes1
        self.modes2 = modes2
        self.modes3 = modes3

        self.scale = 1.0 / (in_channels * out_channels)
        self.weights1 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes1, self.modes2, self.modes3, dtype=torch.cfloat))
        self.weights2 = nn.Parameter(self.scale * torch.rand(in_channels, out_channels, self.modes1, self.modes2, self.modes3, dtype=torch.cfloat))

    def compl_mul3d(self, input, weights):
        return torch.einsum("bixyz,ioxyz->boxyz", input, weights)

    def forward(self, x):
        batchsize = x.shape[0]
        
        # Compute Fourier coefficients
        x_ft = torch.fft.rfftn(x, dim=[-3, -2, -1])

        out_ft = torch.zeros(batchsize, self.out_channels, x.size(-3), x.size(-2), x.size(-1)//2 + 1, dtype=torch.cfloat, device=x.device)
        
        out_ft[:, :, :self.modes1, :self.modes2, :self.modes3] = \
            self.compl_mul3d(x_ft[:, :, :self.modes1, :self.modes2, :self.modes3], self.weights1)
        out_ft[:, :, -self.modes1:, :self.modes2, :self.modes3] = \
            self.compl_mul3d(x_ft[:, :, -self.modes1:, :self.modes2, :self.modes3], self.weights2)

        # Return to spatial domain
        x = torch.fft.irfftn(out_ft, s=(x.size(-3), x.size(-2), x.size(-1)))
        return x


class FNO3d(nn.Module):
    """Full 3D Fourier Neural Operator Network"""
    def __init__(self, in_channels=2, out_channels=1, modes=8, width=20):
        super(FNO3d, self).__init__()
        self.modes = modes
        self.width = width

        self.fc0 = nn.Linear(in_channels, self.width)

        self.conv0 = SpectralConv3d(self.width, self.width, self.modes, self.modes, self.modes)
        self.conv1 = SpectralConv3d(self.width, self.width, self.modes, self.modes, self.modes)
        self.w0 = nn.Conv3d(self.width, self.width, 1)
        self.w1 = nn.Conv3d(self.width, self.width, 1)

        self.fc1 = nn.Linear(self.width, 64)
        self.fc2 = nn.Linear(64, out_channels)

    def forward(self, x):
        # Input shape: (batch, in_channels, dx, dy, dz)
        x = x.permute(0, 2, 3, 4, 1)
        x = self.fc0(x)
        x = x.permute(0, 4, 1, 2, 3)

        x1 = self.conv0(x)
        x2 = self.w0(x)
        x = F.gelu(x1 + x2)

        x1 = self.conv1(x)
        x2 = self.w1(x)
        x = F.gelu(x1 + x2)

        x = x.permute(0, 2, 3, 4, 1)
        x = F.gelu(self.fc1(x))
        x = self.fc2(x)
        x = x.permute(0, 4, 1, 2, 3)
        return x


if __name__ == "__main__":
    model = FNO3d(in_channels=2, out_channels=1, modes=8, width=20)
    dummy_input = torch.randn(2, 2, 32, 32, 32)
    output = model(dummy_input)
    print(f"[Phase 10 FNO Test] Input shape: {dummy_input.shape} -> Output shape: {output.shape}")