import os

import joblib
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
from sklearn.model_selection import GroupShuffleSplit

DATA_PATH = "data/processed/inhA_processed_ml.csv"
METRICS_PATH = "models/scaffold_validation_metrics.csv"
SPLIT_PATH = "models/scaffold_validation_split.csv"
RANDOM_STATE = 42
TEST_FRACTION = 0.20

os.makedirs("models", exist_ok=True)
df = pd.read_csv(DATA_PATH)
feature_cols = [column for column in df.columns if column.startswith("Bit_")]
if not feature_cols:
    raise ValueError(f"No fingerprint columns found in {DATA_PATH}")

X = df[feature_cols].to_numpy()
y = df["class"].to_numpy(dtype=int)
groups = np.asarray(
    [MurckoScaffold.MurckoScaffoldSmiles(smiles=smiles) for smiles in df["canonical_smiles"]]
)

# Select a deterministic scaffold split using only split balance and class presence,
# never model performance, to avoid choosing a favorable test set.
splitter = GroupShuffleSplit(
    n_splits=1000,
    test_size=TEST_FRACTION,
    random_state=RANDOM_STATE,
)
valid_splits = []
for train_idx, test_idx in splitter.split(X, y, groups):
    if len(np.unique(y[train_idx])) < 2 or len(np.unique(y[test_idx])) < 2:
        continue
    test_fraction = len(test_idx) / len(df)
    if not 0.15 <= test_fraction <= 0.30:
        continue
    valid_splits.append((abs(test_fraction - TEST_FRACTION), train_idx, test_idx))

if not valid_splits:
    raise RuntimeError("Could not construct a scaffold-disjoint split containing both classes.")

_, train_idx, test_idx = min(valid_splits, key=lambda split: split[0])
train_scaffolds = set(groups[train_idx])
test_scaffolds = set(groups[test_idx])
if train_scaffolds.intersection(test_scaffolds):
    raise AssertionError("Scaffold leakage detected between train and test partitions.")

model = RandomForestClassifier(
    n_estimators=300,
    max_depth=15,
    class_weight="balanced",
    random_state=RANDOM_STATE,
    n_jobs=-1,
)
model.fit(X[train_idx], y[train_idx])
probabilities = model.predict_proba(X[test_idx])[:, 1]
predictions = (probabilities >= 0.50).astype(int)
tn, fp, fn, tp = confusion_matrix(y[test_idx], predictions, labels=[0, 1]).ravel()

metrics = pd.DataFrame([{
    "Model": "RandomForest",
    "Split_Type": "Bemis-Murcko scaffold holdout",
    "ROC_AUC": roc_auc_score(y[test_idx], probabilities),
    "Accuracy": accuracy_score(y[test_idx], predictions),
    "F1": f1_score(y[test_idx], predictions, zero_division=0),
    "MCC": matthews_corrcoef(y[test_idx], predictions),
    "True_Negatives": tn,
    "False_Positives": fp,
    "False_Negatives": fn,
    "True_Positives": tp,
    "Decision_Threshold": 0.50,
    "Train_Compounds": len(train_idx),
    "Test_Compounds": len(test_idx),
    "Train_Scaffolds": len(train_scaffolds),
    "Test_Scaffolds": len(test_scaffolds),
    "Scaffold_Overlap": len(train_scaffolds.intersection(test_scaffolds)),
    "Split_Random_State": RANDOM_STATE,
    "Candidate_Splits_Checked": 1000,
}])
metrics.to_csv(METRICS_PATH, index=False)

split_output = df[["canonical_smiles", "pIC50", "class"]].copy().reset_index(drop=True)
split_output["murcko_scaffold"] = groups
split_output["partition"] = "train"
split_output.loc[test_idx, "partition"] = "test"
split_output.to_csv(SPLIT_PATH, index=False)

print(metrics.to_string(index=False))
print(f"Metrics saved to {METRICS_PATH}")
print(f"Split assignments saved to {SPLIT_PATH}")
