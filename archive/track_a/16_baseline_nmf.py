import os
import json
import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import NMF
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score

def load_preprocessed_matrices():
    print("[DATA] Initializing data stream from preprocessed arrays...")
    X_raw = np.load("output/nn_X.npy", allow_pickle=True)
    y_raw = np.load("output/nn_y.npy")

    if isinstance(X_raw, np.ndarray) and X_raw.dtype == object:
        X_sparse = sp.vstack(X_raw)
    else:
        X_sparse = sp.csr_matrix(X_raw)

    X_dense = X_sparse.toarray()
    # NMF requires non-negative input; shift log-normalized data to be non-negative
    if X_dense.min() < 0:
        X_dense = X_dense - X_dense.min() + 1e-10
    print(f"[OK] Data Stream Ready. Shape: {X_dense.shape[0]} cells x {X_dense.shape[1]} genes")
    return X_dense, y_raw

def main():
    print("[CPU] Computational Backend: Using CPU (NMF via scikit-learn)")
    print("[PHASE 2] Executing Phase 2: NMF Matrix Factorization & Cellular State Mapping...")
    
    X_dense, y = load_preprocessed_matrices()
    
    print("\n[NMF] Executing NMF Decomposition (n_components=4)...")
    nmf_model = NMF(
        n_components=4,
        init='nndsvda',
        max_iter=500,
        random_state=42,
        solver='cd',
        beta_loss='frobenius'
    )
    W = nmf_model.fit_transform(X_dense)
    H = nmf_model.components_.T
    
    print(f"[OK] NMF Decomposition Complete.")
    print(f"   W (Cell x State) shape: {W.shape}")
    print(f"   H (Gene x State) shape: {H.shape}")
    print(f"   Reconstruction Error (Frobenius): {nmf_model.reconstruction_err_:.4f}")
    print(f"   Meta-Modules: AC (0), MES (1), NPC (2), OPC (3)")

    print("\n[PHASE 2B] Executing Phase 2B: Frozen Fraction Linear Evaluation Head...")
    X_train, X_test, y_train, y_test = train_test_split(
        W, y, test_size=0.20, random_state=42, stratify=y
    )

    clf = LogisticRegression(
        max_iter=1000,
        random_state=42,
        solver='lbfgs',
        C=1.0
    )
    clf.fit(X_train, y_train)

    print("\n[EVAL] Evaluating NMF Baseline on Test Partition...")
    test_predictions = clf.predict(X_test)
    test_probabilities = clf.predict_proba(X_test)

    accuracy = accuracy_score(y_test, test_predictions)
    f1 = f1_score(y_test, test_predictions, average='macro')
    auc = roc_auc_score(y_test, test_probabilities, multi_class='ovr', average='macro')

    print("\n================== [NMF] NMF MATRIX FACTORIZATION BASELINE RESULTS ==================")
    print(f"Testing Classification Accuracy: {accuracy * 100:.2f}%")
    print(f"Macro-Averaged F1 Score:         {f1:.4f}")
    print(f"Macro-Averaged AUC (OvR):        {auc:.4f}")
    print("=========================================================================")

    os.makedirs("output", exist_ok=True)
    metrics_export = {"accuracy": accuracy, "macro_f1": f1, "macro_auc": auc}
    with open("output/nmf_metrics.json", "w") as f:
        json.dump(metrics_export, f, indent=4)

    np.save("output/nmf_fractions.npy", W)
    print("[OK] Performance metrics saved to output/nmf_metrics.json")
    print("[OK] Cell-state fraction matrix saved to output/nmf_fractions.npy")

if __name__ == "__main__":
    main()