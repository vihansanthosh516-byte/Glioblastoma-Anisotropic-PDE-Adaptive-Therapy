import os
import json
import numpy as np
import scipy.sparse as sp
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score

# =====================================================================
# 1. HARDWARE ACCELERATION SETUP
# =====================================================================
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Computational Backend: Using {DEVICE.type.upper()}")
if DEVICE.type == 'cuda':
    print(f"   Target Device: {torch.cuda.get_device_name(0)}")
    torch.backends.cudnn.benchmark = True

# =====================================================================
# 2. DATA PROCESSING & STORAGE MANAGEMENT
# =====================================================================
class SingleCellDataset(Dataset):
    def __init__(self, X_sparse, y):
        # Keep data in CPU RAM as sparse CSR format to prevent RAM explosion
        self.X = X_sparse.tocsr()
        self.y = torch.LongTensor(y)
        
    def __len__(self):
        return self.X.shape[0]
        
    def __getitem__(self, idx):
        # Dynamically realize a dense row matrix right before shipping to GPU
        x_dense = torch.FloatTensor(self.X[idx].toarray().squeeze())
        return x_dense, self.y[idx]

def load_preprocessed_matrices():
    print("Initializing data stream from preprocessed arrays...")
    X_raw = np.load("output/nn_X.npy", allow_pickle=True)
    y_raw = np.load("output/nn_y.npy")
    
    if isinstance(X_raw, np.ndarray) and X_raw.dtype == object:
        X_sparse = sp.vstack(X_raw)
    else:
        X_sparse = sp.csr_matrix(X_raw)
        
    print(f"Data Stream Ready. Shape: {X_sparse.shape[0]} cells x {X_sparse.shape[1]} genes")
    return X_sparse, y_raw

# =====================================================================
# 3. CORE ARCHITECTURE: SPECIFIC TRANSCRIPTOMIC VAE
# =====================================================================
class SingleCellVAE(nn.Module):
    def __init__(self, input_dim, latent_dim=32):
        super(SingleCellVAE, self).__init__()
        
        # Encoder Module (2,500 HVGs -> 256 -> 128 -> Distributions)
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.ReLU()
        )
        self.fc_mu = nn.Linear(128, latent_dim)
        self.fc_logvar = nn.Linear(128, latent_dim)
        
        # Decoder Module (32 -> 128 -> 256 -> 2,500 HVG Reconstructions)
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Linear(128, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Linear(256, input_dim)
        )
        
    def encode(self, x):
        h = self.encoder(x)
        return self.fc_mu(h), self.fc_logvar(h)
        
    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std
        
    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        return self.decoder(z), mu, logvar

# Mathematical Optimization: Evidence Lower Bound (ELBO) Loss Function
def compute_elbo_loss(recon_x, x, mu, logvar):
    # Reconstruction Loss modeled via Mean Squared Error on regularized profiles
    recon_loss = nn.functional.mse_loss(recon_x, x, reduction='sum')
    # Analytical Kullback-Leibler Divergence for Normal Gaussian Prior
    kl_divergence = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    return recon_loss + kl_divergence

# =====================================================================
# 4. FROZEN-LATENT LINEAR EVALUATION HEAD
# =====================================================================
class LinearEvaluationHead(nn.Module):
    def __init__(self, latent_dim=32, num_classes=3):
        super(LinearEvaluationHead, self).__init__()
        self.classifier = nn.Linear(latent_dim, num_classes)
        
    def forward(self, z):
        return self.classifier(z)

# =====================================================================
# 5. PIPELINE EXECUTION
# =====================================================================
def main():
    X_sparse, y = load_preprocessed_matrices()
    input_dim = X_sparse.shape[1]
    
    # Stratified 80/20 train/test distribution split
    indices = np.arange(X_sparse.shape[0])
    train_idx, test_idx, y_train, y_test = train_test_split(
        indices, y, test_size=0.20, random_state=42, stratify=y
    )
    
    train_loader = DataLoader(SingleCellDataset(X_sparse[train_idx], y_train), batch_size=256, shuffle=True)
    test_loader = DataLoader(SingleCellDataset(X_sparse[test_idx], y_test), batch_size=256, shuffle=False)
    
    # Initialize and target GPU acceleration
    vae = SingleCellVAE(input_dim=input_dim, latent_dim=32).to(DEVICE)
    evaluation_head = LinearEvaluationHead(latent_dim=32, num_classes=3).to(DEVICE)
    
    vae_optimizer = optim.Adam(vae.parameters(), lr=0.001)
    clf_optimizer = optim.Adam(evaluation_head.parameters(), lr=0.005)
    classification_criterion = nn.CrossEntropyLoss()
    
    # --- PHASE 1A: UNSUPERVISED VAE OPTIMIZATION ---
    print("\nExecuting Phase 1A: Unsupervised VAE Latent Compression...")
    vae.train()
    for epoch in range(20):
        running_elbo = 0
        for batch_x, _ in train_loader:
            batch_x = batch_x.to(DEVICE)
            vae_optimizer.zero_grad()
            recon_x, mu, logvar = vae(batch_x)
            loss = compute_elbo_loss(recon_x, batch_x, mu, logvar)
            loss.backward()
            vae_optimizer.step()
            running_elbo += loss.item()
        print(f"   Epoch {epoch+1:02d}/20 | Normalized ELBO Loss: {running_elbo / len(train_loader.dataset):.4f}")
        
    # --- PHASE 1B: EXTRACT FROZEN EMBEDDINGS & TRAIN EVAL HEAD ---
    print("\nFreezing VAE Latent Parameters & Optimizing Evaluation Head...")
    vae.eval()
    
    for epoch in range(15):
        evaluation_head.train()
        running_clf_loss = 0
        for batch_x, batch_y in train_loader:
            batch_x, batch_y = batch_x.to(DEVICE), batch_y.to(DEVICE)
            clf_optimizer.zero_grad()
            with torch.no_grad():
                mu, logvar = vae.encode(batch_x) # Utilizing structural mean space
            outputs = evaluation_head(mu)
            loss = classification_criterion(outputs, batch_y)
            loss.backward()
            clf_optimizer.step()
            running_clf_loss += loss.item()
            
    # --- PHASE 1C: GRADIENT EVALUATION MATRIX ---
    print("\nEvaluating Deep Generative scVI Baseline on Test Partition...")
    evaluation_head.eval()
    test_predictions = []
    test_probabilities = []
    
    with torch.no_grad():
        for batch_x, _ in test_loader:
            batch_x = batch_x.to(DEVICE)
            mu, logvar = vae.encode(batch_x)
            outputs = evaluation_head(mu)
            probs = nn.functional.softmax(outputs, dim=1)
            preds = torch.argmax(probs, dim=1)
            test_predictions.extend(preds.cpu().numpy())
            test_probabilities.extend(probs.cpu().numpy())
            
    test_predictions = np.array(test_predictions)
    test_probabilities = np.array(test_probabilities)
    
    # Derive Formal Computational Metrics
    accuracy = accuracy_score(y_test, test_predictions)
    f1 = f1_score(y_test, test_predictions, average='macro')
    auc = roc_auc_score(y_test, test_probabilities, multi_class='ovr', average='macro')
    
    print("\n================== scVI GENE MATRIX BASELINE RESULTS ==================")
    print(f"Testing Classification Accuracy: {accuracy * 100:.2f}%")
    print(f"Macro-Averaged F1 Score:         {f1:.4f}")
    print(f"Macro-Averaged AUC (OvR):        {auc:.4f}")
    print("=========================================================================")
    
    # Export metrics metadata
    os.makedirs("output", exist_ok=True)
    metrics_export = {"accuracy": accuracy, "macro_f1": f1, "macro_auc": auc}
    with open("output/scvi_metrics.json", "w") as f:
        json.dump(metrics_export, f, indent=4)
        
    # Extract and persist full dataset latent embeddings for Phase 3 diagnostics
    print("\nExtracting and compiling global 140k latent coordinates for Phase 3...")
    global_loader = DataLoader(SingleCellDataset(X_sparse, y), batch_size=512, shuffle=False)
    global_latent_matrix = []
    
    with torch.no_grad():
        for batch_x, _ in global_loader:
            batch_x = batch_x.to(DEVICE)
            mu, _ = vae.encode(batch_x)
            global_latent_matrix.extend(mu.cpu().numpy())
            
    np.save("output/scvi_latent.npy", np.array(global_latent_matrix))
    print("Complete latent matrix successfully saved at output/scvi_latent.npy!")

if __name__ == "__main__":
    main()