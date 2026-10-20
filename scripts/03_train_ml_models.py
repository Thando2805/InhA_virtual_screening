import os
import joblib
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, matthews_corrcoef, confusion_matrix

try:
    from xgboost import XGBClassifier
except ImportError:
    XGBClassifier = None

os.makedirs("models", exist_ok=True)

# 1. Load dataset
df = pd.read_csv("data/processed/inhA_processed_ml.csv")
feature_cols = [col for col in df.columns if col.startswith("Bit_")]
X = df[feature_cols].values
y = df['class'].values

num_pos = sum(y == 1)
num_neg = sum(y == 0)
scale_pos_weight_value = num_neg / num_pos if num_pos > 0 else 1.0

print(f"Loaded {len(df)} compounds. Actives: {num_pos} | Inactives: {num_neg}")
print(f"Applying Class Weighting Ratio: {scale_pos_weight_value:.2f}")

# 2. Train/Test Split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.20, random_state=42, stratify=y
)

# 3. Classifiers with Balanced Weights
models = {
    "RandomForest": RandomForestClassifier(
        n_estimators=300, 
        max_depth=15, 
        class_weight="balanced", # Penalizes missing active compounds
        random_state=42, 
        n_jobs=-1
    ),
}

if XGBClassifier is not None:
    models["XGBoost"] = XGBClassifier(
        n_estimators=300, 
        max_depth=6, 
        learning_rate=0.03, 
        scale_pos_weight=scale_pos_weight_value, # Adjusts for class imbalance
        random_state=42, 
        eval_metric="logloss"
    )

best_model = None
best_roc_auc = -np.inf
best_model_name = ""
metrics = []

for name, model in models.items():
    model.fit(X_train, y_train)
    y_proba = model.predict_proba(X_test)[:, 1]
    y_pred = (y_proba >= 0.50).astype(int)

    roc = roc_auc_score(y_test, y_proba)
    accuracy = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, zero_division=0)
    mcc = matthews_corrcoef(y_test, y_pred)
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred, labels=[0, 1]).ravel()

    print(f"\n--- {name} (Balanced) ---")
    print(f"ROC-AUC: {roc:.4f} | Accuracy: {accuracy:.4f} | F1: {f1:.4f} | MCC: {mcc:.4f}")
    metrics.append({
        "Model": name,
        "ROC_AUC": roc,
        "Accuracy": accuracy,
        "F1": f1,
        "MCC": mcc,
        "True_Negatives": tn,
        "False_Positives": fp,
        "False_Negatives": fn,
        "True_Positives": tp,
        "Decision_Threshold": 0.50,
        "Train_Compounds": len(X_train),
        "Test_Compounds": len(X_test),
        "Split_Random_State": 42,
    })

    if roc > best_roc_auc:
        best_roc_auc = roc
        best_model = model
        best_model_name = name

# Save updated model
joblib.dump(best_model, "models/best_inhA_model.joblib")
print(f"\nRetrained {best_model_name} exported to models/best_inhA_model.joblib")

metrics_path = "models/training_metrics.csv"
pd.DataFrame(metrics).to_csv(metrics_path, index=False)
print(f"Held-out evaluation metrics saved to {metrics_path}")
if XGBClassifier is None:
    print("XGBoost was unavailable; only RandomForest was trained and evaluated.")