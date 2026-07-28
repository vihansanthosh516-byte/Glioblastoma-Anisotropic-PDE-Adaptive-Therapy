#!/usr/bin/env python3
"""
Phase 10: FNO Training Pipeline
================================
Train Fourier Neural Operator to accelerate 3D tumor growth PDE.
"""
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from pathlib import Path
import sys
import json

sys.path.insert(0, str(Path(__file__).parent))
from fno_solver import FNO3d


def train_fno(
    data_path="output/fno_dataset.pt",
    save_path="output/fno_model.pth",
    epochs=50,
    batch_size=4,
    lr=1e-3,
    modes=8,
    width=20,
    device=None
):
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[Phase 10 FNO Training] Using device: {device}")

    # Load dataset
    data = torch.load(data_path, map_location=device)
    x_data = data["x"]
    y_data = data["y"]
    print(f"[Data] X: {x_data.shape}, Y: {y_data.shape}")

    # Split train/val
    n = len(x_data)
    train_n = int(0.8 * n)
    train_dataset = TensorDataset(x_data[:train_n], y_data[:train_n])
    val_dataset = TensorDataset(x_data[train_n:], y_data[train_n:])

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

    # Model
    model = FNO3d(in_channels=2, out_channels=1, modes=8, width=20).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=5, factor=0.5)
    criterion = nn.MSELoss()

    print(f"[Training] Epochs: {epochs}, Batch: {batch_size}, LR: {lr}")
    print(f"[Training] Train: {len(train_dataset)}, Val: {len(val_dataset)}")

    best_val_loss = float('inf')
    history = {"train_loss": [], "val_loss": []}

    for epoch in range(1, epochs + 1):
        model.train()
        train_loss = 0.0
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            pred = model(xb)
            loss = criterion(pred, yb)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * xb.size(0)
        train_loss /= len(train_dataset)

        # Validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for xb, yb in val_loader:
                xb, yb = xb.to(device), yb.to(device)
                pred = model(xb)
                loss = criterion(pred, yb)
                val_loss += loss.item() * xb.size(0)
        val_loss /= len(val_dataset)

        scheduler.step(val_loss)
        history["train_loss"].append(train_loss)
        history["val_loss"].append(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), save_path)
            status = " *** BEST ***"
        else:
            status = ""

        if epoch % 5 == 0 or epoch == 1 or epoch == epochs:
            print(f"  Epoch {epoch:3d}/{epochs}: Train MSE={train_loss:.6f}, Val MSE={val_loss:.6f}{status}")

    # Save history
    torch.save({"model_state": model.state_dict(), "history": history, "config": {"modes": 8, "width": 20}}, save_path)
    print(f"[Training] Best val MSE: {best_val_loss:.6f}. Model saved to {save_path}")
    return model, history


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Train FNO for tumor growth PDE")
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--data", type=str, default="output/fno_dataset.pt")
    parser.add_argument("--save", type=str, default="output/fno_model.pth")
    args = parser.parse_args()

    train_fno(
        data_path=args.data,
        save_path=args.save,
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr
    )