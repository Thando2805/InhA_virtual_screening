import os

import numpy as np
import pandas as pd
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold

DATA_PATH = "data/processed/inhA_processed_ml.csv"
FOLDS_PATH = "models/repeated_scaffold_cv_folds.csv"
SUMMARY_PATH = "models/repeated_scaffold_cv_summary.csv"
N_SPLITS = 5
N_REPEATS = 5
BASE_RANDOM_STATE = 42
THRESHOLD = 0.50

os.makedirs("models", exist_ok=True)
df = pd.read_csv(DATA_PATH)
feature_cols = [column for column in df.columns if column.startswith("Bit_")]
if not feature_cols:
    raise ValueError(f"No fingerprint columns found in {DATA_PATH}")

X = df[feature_cols].to_numpy()
y = df["class"].to_numpy(dtype=int)
groups = np.asarray([
    MurckoScaffold.MurckoScaffoldSmiles(smiles=smiles)
    for smiles in df["canonical_smiles"]
])

fold_rows = []
for repeat in range(N_REPEATS):
    random_state = BASE_RANDOM_STATE + repeat
    splitter = StratifiedGroupKFold(
        n_splits=N_SPLITS,
        shuffle=True,
        random_state=random_state,
    )

    for fold, (train_idx, test_idx) in enumerate(splitter.split(X, y, groups), start=1):
        train_scaffolds = set(groups[train_idx])
        test_scaffolds = set(groups[test_idx])
        scaffold_overlap = train_scaffolds.intersection(test_scaffolds)
        if scaffold_overlap:
            raise AssertionError(
                f"Scaffold leakage in repeat {repeat + 1}, fold {fold}: "
                f"{len(scaffold_overlap)} shared scaffolds"
            )
        if len(np.unique(y[test_idx])) < 2:
            raise RuntimeError(
                f"Test fold {fold} in repeat {repeat + 1} contains only one class; "
                "cannot calculate ROC-AUC."
            )

        model = RandomForestClassifier(
            n_estimators=300,
            max_depth=15,
            class_weight="balanced",
            random_state=random_state * N_SPLITS + fold,
            n_jobs=-1,
        )
        model.fit(X[train_idx], y[train_idx])
        probabilities = model.predict_proba(X[test_idx])[:, 1]
        predictions = (probabilities >= THRESHOLD).astype(int)
        tn, fp, fn, tp = confusion_matrix(y[test_idx], predictions, labels=[0, 1]).ravel()

        fold_rows.append({
            "Repeat": repeat + 1,
            "Fold": fold,
            "Split_Random_State": random_state,
            "ROC_AUC": roc_auc_score(y[test_idx], probabilities),
            "Accuracy": accuracy_score(y[test_idx], predictions),
            "F1": f1_score(y[test_idx], predictions, zero_division=0),
            "MCC": matthews_corrcoef(y[test_idx], predictions),
            "True_Negatives": tn,
            "False_Positives": fp,
            "False_Negatives": fn,
            "True_Positives": tp,
            "Decision_Threshold": THRESHOLD,
            "Train_Compounds": len(train_idx),
            "Test_Compounds": len(test_idx),
            "Train_Scaffolds": len(train_scaffolds),
            "Test_Scaffolds": len(test_scaffolds),
            "Scaffold_Overlap": len(scaffold_overlap),
        })

fold_metrics = pd.DataFrame(fold_rows)
fold_metrics.to_csv(FOLDS_PATH, index=False)

metric_columns = ["ROC_AUC", "Accuracy", "F1", "MCC"]
summary_rows = []
for metric in metric_columns:
    values = fold_metrics[metric]
    summary_rows.append({
        "Metric": metric,
        "Mean": values.mean(),
        "SD": values.std(ddof=1),
        "Min": values.min(),
        "Max": values.max(),
        "Fold_Count": len(values),
        "Repeats": N_REPEATS,
        "Folds_Per_Repeat": N_SPLITS,
        "Validation": "Repeated 5-fold Bemis-Murcko StratifiedGroupKFold",
        "Base_Random_State": BASE_RANDOM_STATE,
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(SUMMARY_PATH, index=False)

if not (fold_metrics["Scaffold_Overlap"] == 0).all():
    raise AssertionError("One or more saved folds have scaffold overlap.")

print("Repeated 5-fold scaffold-aware cross-validation (mean +/- SD):")
for row in summary.itertuples(index=False):
    print(f"{row.Metric}: {row.Mean:.4f} +/- {row.SD:.4f} "
          f"(range {row.Min:.4f}-{row.Max:.4f})")
print(f"Folds saved to {FOLDS_PATH}")
print(f"Summary saved to {SUMMARY_PATH}")
print(f"Verified scaffold overlap across {len(fold_metrics)} folds: 0")
