import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Ensure output directory exists
os.makedirs("data/processed", exist_ok=True)

# 1. Load raw ChEMBL activity data
df = pd.read_csv("data/raw/inhA_chembl_ic50_raw.csv")
print(f"Loaded {len(df)} raw entries.")

# Keep only directly comparable, uncensored IC50 measurements in nM.
required_columns = {
    'standard_type', 'standard_units', 'standard_relation',
    'standard_value', 'canonical_smiles'
}
missing_columns = required_columns.difference(df.columns)
if missing_columns:
    raise ValueError(f"Raw activity file is missing required columns: {sorted(missing_columns)}")

raw_count = len(df)
df['standard_type'] = df['standard_type'].astype('string').str.strip().str.upper()
df['standard_units'] = df['standard_units'].astype('string').str.strip().str.lower()
df['standard_relation'] = df['standard_relation'].astype('string').str.strip()
df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
df = df[
    (df['standard_type'] == 'IC50')
    & (df['standard_units'] == 'nm')
    & (df['standard_relation'] == '=')
    & (df['standard_value'] > 0)
].dropna(subset=['standard_value', 'canonical_smiles'])
eligible_count = len(df)

# Convert exact IC50 values from nM to molar concentration before taking -log10.
df['pIC50'] = -np.log10(df['standard_value'] * 1e-9)

# Aggregate replicate measurements first, then apply the 1 uM threshold.
df = df.groupby('canonical_smiles', as_index=False).agg({'pIC50': 'median'})
df['class'] = (df['pIC50'] >= 6.0).astype(int)
print(f"Dataset size after deduplication: {len(df)} unique compounds")
print(f"Class distribution:\n{df['class'].value_counts()}")

# 5. Generate 2048-bit Morgan Fingerprints (ECFP4)
def smiles_to_ecfp4(smiles, radius=2, nBits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    return list(fp)

print("Computing ECFP4 fingerprints...")
fingerprints = []
valid_indices = []

for idx, smiles in enumerate(df['canonical_smiles']):
    fp = smiles_to_ecfp4(smiles)
    if fp is not None:
        fingerprints.append(fp)
        valid_indices.append(idx)

# Filter dataframe to valid RDKit structures
df_clean = df.iloc[valid_indices].copy()
fp_cols = [f"Bit_{i}" for i in range(2048)]
fp_df = pd.DataFrame(fingerprints, columns=fp_cols, index=df_clean.index)

# Combine SMILES, target metrics, and fingerprint feature matrix
processed_df = pd.concat([df_clean[['canonical_smiles', 'pIC50', 'class']], fp_df], axis=1)

# Export processed dataset
output_path = "data/processed/inhA_processed_ml.csv"
processed_df.to_csv(output_path, index=False)
print(f"Successfully exported ML-ready dataset to {output_path}")

audit = pd.DataFrame([
    {'Metric': 'raw_activity_records', 'Count': raw_count},
    {'Metric': 'eligible_exact_ic50_nM_records', 'Count': eligible_count},
    {'Metric': 'excluded_records', 'Count': raw_count - eligible_count},
    {'Metric': 'unique_valid_structures', 'Count': len(processed_df)},
    {'Metric': 'active_structures', 'Count': int((processed_df['class'] == 1).sum())},
    {'Metric': 'inactive_structures', 'Count': int((processed_df['class'] == 0).sum())},
])
audit.to_csv('data/processed/inhA_dataset_audit.csv', index=False)
print("Dataset filtering audit saved to data/processed/inhA_dataset_audit.csv")