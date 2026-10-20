import os
import joblib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Complete Zimbabwean Phytochemical Candidate Library (21 Proposal Plants + Lantana + Control)
candidates_data = [
    # --- Native / Indigenous Species (Part A) ---
    {"Compound_Name": "Carissone", "Plant_Source": "Carissa edulis / C. bispinosa", "SMILES": "CC1=CC(=O)[C@@H]2C[C@@H](C(C)(C)O)[C@H](C)[C@@]2(C)C1"},
    {"Compound_Name": "Ursolic Acid", "Plant_Source": "Syzygium cordatum", "SMILES": "CC1CCC2(CCC3(C(=CCC4C3(C)CCC5C4(C)CCC(O)C5(C)C)C2C1C)C)C(=O)O"},
    {"Compound_Name": "Lupeol", "Plant_Source": "Ficus sycomorus", "SMILES": "CC(=C)[C@H]1CC[C@@]2(C)[C@@H]1CC[C@]3(C)[C@@H]2CC=C4[C@@H]3C(C)(C)[C@@H](O)CC4"},
    {"Compound_Name": "Beta-Sitosterol", "Plant_Source": "Coleochloa setifera", "SMILES": "CCC(CC(C)C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C"},
    {"Compound_Name": "Securinine", "Plant_Source": "Flueggea virosa", "SMILES": "O=C1C=C2[C@@H]3N4CCC[C@H]4[C@@H](C1)C2=CO3"},
    {"Compound_Name": "Sarsasapogenin", "Plant_Source": "Asparagus africanus", "SMILES": "CC1CCC2C(C1)C3C(O2)CC4C3(CCC5C4CC=C6CC(O)CCC56C)C"},
    {"Compound_Name": "Stigmasterol", "Plant_Source": "Annona stenophylla", "SMILES": "CCC(C=CC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C)C(C)C"},
    {"Compound_Name": "Kaempferol", "Plant_Source": "Pellaea sp.", "SMILES": "C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O"},
    {"Compound_Name": "Verbascoside", "Plant_Source": "Lippia javanica", "SMILES": "O=C(/C=C/c1ccc(O)c(O)c1)O[C@@H]2[C@@H](O)[C@H](O[C@@H](O[C@H]3[C@@H](O)[C@H](O)[C@@H](O)O3)C2)OCCOc4ccc(O)c(O)c4"},
    {"Compound_Name": "Piperitenone", "Plant_Source": "Lippia javanica", "SMILES": "CC1=CC(=O)C(=C(C)C)CC1"},
    {"Compound_Name": "Camphor", "Plant_Source": "Lippia javanica", "SMILES": "CC1(C)C2CCC1(C)C(=O)C2"},
    {"Compound_Name": "Quercetin", "Plant_Source": "Myrothamnus flabellifolius", "SMILES": "O=C1C(O)=C(c2ccc(O)c(O)c2)Oc3cc(O)cc(O)c13"},
    {"Compound_Name": "Arbutin", "Plant_Source": "Myrothamnus flabellifolius", "SMILES": "C1=CC(=CC=C1O)OC2C(C(C(C(O2)CO)O)O)O"},
    {"Compound_Name": "Quercitrin", "Plant_Source": "Piliostigma thonningii", "SMILES": "CC1C(C(C(C(O1)OC2=C(OC3=CC(=CC(=C3C2=O)O)O)C4=CC(=C(C=C4)O)O)O)O)O"},
    {"Compound_Name": "Neorautenanol", "Plant_Source": "Dalbergia melanoxylon", "SMILES": "CC1(C)C=CC2=C(O1)C=CC3=C2C(C4=C(O3)C=CC(=C4)O)O"},
    {"Compound_Name": "Amentoflavone", "Plant_Source": "Searsia pyroides", "SMILES": "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)C4=C(C=CC(=C4)C5=CC(=O)C6=C(O5)C=C(C=C6O)O)O)O"},
    {"Compound_Name": "7-Methyljuglone", "Plant_Source": "Euclea crispa", "SMILES": "CC1=CC2=C(C(=O)C=CC2=O)C(=C1)O"},
    {"Compound_Name": "Friedelin", "Plant_Source": "Bridelia mollis", "SMILES": "CC1CC2(CCC3(C(=CCC4C3(C)CCC5C4(C)CCC(=O)C5(C)C)C2C1C)C)C"},
    {"Compound_Name": "Terpineol", "Plant_Source": "Vitex payos", "SMILES": "CC1=CCC(CC1)C(C)(C)O"},
    {"Compound_Name": "Rutin", "Plant_Source": "Ficus ingens", "SMILES": "CC1C(C(C(C(O1)OCC2C(C(C(C(O2)OC3=C(OC4=CC(=CC(=C4C3=O)O)O)C5=CC(=C(C=C5)O)O)O)O)O)O)O)O"},
    
    # --- Formulated / Naturalized Species ---
    {"Compound_Name": "Lantadene A", "Plant_Source": "Lantana camara", "SMILES": "CC(=C/C(=O)O[C@H]1CC2(C)C(=CCC3C2(C)CCC4C3(C)CCC(C)(C(=O)O)C4)C1(C)C)/C"},
    {"Compound_Name": "Guaijaverin", "Plant_Source": "Psidium guajava", "SMILES": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)OC4C(C(C(O4)O)O)O)O)O"},
    {"Compound_Name": "Eucalyptol (1,8-Cineole)", "Plant_Source": "Eucalyptus camaldulensis", "SMILES": "CC1(C2CCC(O1)(CC2)C)C"},
    {"Compound_Name": "Limonene", "Plant_Source": "Citrus limon", "SMILES": "CC1=CCC(CC1)C(=C)C"},

    # --- Positive Control ---
    {"Compound_Name": "PT70 (Control Inhibitor)", "Plant_Source": "Synthetic Control", "SMILES": "CCCCCC1=CC(=C(C=C1)OC2=CC=CC=C2C)O"}
]

df_candidates = pd.DataFrame(candidates_data)

# 1. Generate Fingerprints
def smiles_to_ecfp4(smiles, radius=2, nBits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return list(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits))

fingerprints = [smiles_to_ecfp4(s) for s in df_candidates['SMILES']]
fp_cols = [f"Bit_{i}" for i in range(2048)]
X_screen = pd.DataFrame(fingerprints, columns=fp_cols)

# 2. Screen through trained Machine Learning model
model = joblib.load("models/best_inhA_model.joblib")
probabilities = model.predict_proba(X_screen.values)[:, 1]

df_candidates['Inhibition_Probability'] = probabilities
df_candidates['ML_Prediction'] = (df_candidates['Inhibition_Probability'] >= 0.50).astype(int)

# 3. Rank and Export
df_ranked = df_candidates.sort_values(by="Inhibition_Probability", ascending=False).reset_index(drop=True)

print(f"\n================ Screened {len(df_ranked)} Total Phytochemical Candidates ================")
for idx, row in df_ranked.iterrows():
    status = "ACTIVE HIT" if row['ML_Prediction'] == 1 else "INACTIVE"
    print(f"Rank {idx+1:02d}: {row['Compound_Name']:<20} ({row['Plant_Source']:<30}) | Prob: {row['Inhibition_Probability']:.4f} | [{status}]")

output_path = "data/processed/zimbabwean_phytochemical_hits.csv"
df_ranked.to_csv(output_path, index=False)
print(f"\nFull ranked results saved to {output_path}")