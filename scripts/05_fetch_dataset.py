import os
import time
import requests
import pandas as pd
from rdkit import Chem

# ==============================================================================
# 1. BASELINE SEED DATASET (Covers all 21 species + PT70 control)
# ==============================================================================
SEED_DATA = [
    # 1. Carissa edulis
    {"Compound_Name": "Carissone", "Plant_Source": "Carissa edulis", "SMILES": "CC1=CC(=O)[C@@H]2C[C@@H](C(C)(C)O)[C@H](C)[C@@]2(C)C1"},
    {"Compound_Name": "Eudesmol", "Plant_Source": "Carissa edulis", "SMILES": "CC1=CCCC2(C)C1CCC2(C)C(C)(C)O"},
    
    # 2. Syzygium cordatum
    {"Compound_Name": "Ursolic Acid", "Plant_Source": "Syzygium cordatum", "SMILES": "CC1CCC2(CCC3(C(=CCC4C3(C)CCC5C4(C)CCC(O)C5(C)C)C2C1C)C)C(=O)O"},
    {"Compound_Name": "Betulinic Acid", "Plant_Source": "Syzygium cordatum", "SMILES": "CC(=C)C1CCC2(C1C3CCC4C5(CCC(C(C5CCC4(C3(CC2)C)C)(C)C)O)C)C(=O)O"},
    
    # 3. Ficus sycomorus
    {"Compound_Name": "Sycomoric Acid", "Plant_Source": "Ficus sycomorus", "SMILES": "CC1(C)C2CCC3(C)C(=CCC4C3(C)CCC5C4(C)CCC(O)C5(C)C)C2CC1=O"},
    {"Compound_Name": "Psoralen", "Plant_Source": "Ficus sycomorus", "SMILES": "O=C1C=CC2=CC3=C(C=C2O1)OC=C3"},
    
    # 4. Carissa bispinosa
    {"Compound_Name": "Carisol", "Plant_Source": "Carissa bispinosa", "SMILES": "CC1=CC(=O)C2CC(C(C)(C)O)C(C)C2(C)C1"},
    {"Compound_Name": "Carissanol", "Plant_Source": "Carissa bispinosa", "SMILES": "CC1(C)CCCC2(C)C1CCC3(C)C2CC=C4C3(C)CCC(O)C4(C)C"},
    
    # 5. Coleochloa setifera
    {"Compound_Name": "Luteolin", "Plant_Source": "Coleochloa setifera", "SMILES": "C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O"},
    {"Compound_Name": "Apigenin", "Plant_Source": "Coleochloa setifera", "SMILES": "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O"},
    
    # 6. Flueggea virosa
    {"Compound_Name": "Securinine", "Plant_Source": "Flueggea virosa", "SMILES": "O=C1C=C2[C@@H]3N4CCC[C@H]4[C@@H](C1)C2=CO3"},
    {"Compound_Name": "Virosine A", "Plant_Source": "Flueggea virosa", "SMILES": "CC1=CC2=C(C=C1)N3CCC4(C(=O)OCC34)C2"},
    
    # 7. Asparagus africanus
    {"Compound_Name": "Asparagoside A", "Plant_Source": "Asparagus africanus", "SMILES": "CC1CCC2(C(O1)C(C3C2C4C(=CC5C4(CCC(C5)O)C)CC3)O)C"},
    {"Compound_Name": "Sarsasapogenin", "Plant_Source": "Asparagus africanus", "SMILES": "CC1CCC2C(C1)C3C(O2)CC4C3(CCC5C4CCC6C5(CCC(C6)O)C)C"},
    
    # 8. Annona stenophylla
    {"Compound_Name": "Stenophyllin A", "Plant_Source": "Annona stenophylla", "SMILES": "CCCCCCCCCCCC1CC(O1)CC2CC(O2)CC3CC(=O)OC3"},
    {"Compound_Name": "Anonaine", "Plant_Source": "Annona stenophylla", "SMILES": "C1CNCC2=C1C3=C(C=C2)OCO3"},
    
    # 9. Pellaea sp.
    {"Compound_Name": "Kaempferol", "Plant_Source": "Pellaea sp.", "SMILES": "C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O"},
    {"Compound_Name": "Quercitrin", "Plant_Source": "Pellaea sp.", "SMILES": "CC1C(C(C(C(O1)OC2=C(OC3=CC(=CC(=C3C2=O)O)O)C4=CC(=C(C=C4)O)O)O)O)O"},
    
    # 10. Lippia javanica
    {"Compound_Name": "Verbascoside", "Plant_Source": "Lippia javanica", "SMILES": "CC1C(C(C(C(O1)OCC2C(C(C(C(O2)OCC3=CC(=C(C=C3)O)O)O)OC(=O)/C=C/C4=CC(=C(C=C4)O)O)O)O)O)O"},
    {"Compound_Name": "Piperitenone", "Plant_Source": "Lippia javanica", "SMILES": "CC1=CC(=O)C(=C(C)C)CC1"},
    {"Compound_Name": "Camphor", "Plant_Source": "Lippia javanica", "SMILES": "CC1(C)C2CCC1(C)C(=O)C2"},
    
    # 11. Myrothamnus flabellifolius
    {"Compound_Name": "3,4,5-Tri-O-galloylquinic Acid", "Plant_Source": "Myrothamnus flabellifolius", "SMILES": "C1C(C(C(CC1(C(=O)O)O)OC(=O)C2=CC(=C(C(=C2)O)O)O)OC(=O)C3=CC(=C(C(=C3)O)O)O)OC(=O)C4=CC(=C(C(=C4)O)O)O"},
    {"Compound_Name": "Arbutin", "Plant_Source": "Myrothamnus flabellifolius", "SMILES": "C1=CC(=CC=C1O)OC2C(C(C(C(O2)CO)O)O)O"},
    
    # 12. Piliostigma thonningii
    {"Compound_Name": "Piliostigmin", "Plant_Source": "Piliostigma thonningii", "SMILES": "CC(=O)OC1C(C2=C(C=C1)C(=O)C3=C(O2)C=CC(=C3O)O)O"},
    {"Compound_Name": "Quercetin", "Plant_Source": "Piliostigma thonningii", "SMILES": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O"},
    
    # 13. Dalbergia melanoxylon
    {"Compound_Name": "Melanoxanthin", "Plant_Source": "Dalbergia melanoxylon", "SMILES": "COC1=C(C=C(C=C1)C2=C(C(=O)C3=C(O2)C=C(C=C3O)O)O)O"},
    {"Compound_Name": "Dalbergin", "Plant_Source": "Dalbergia melanoxylon", "SMILES": "COC1=C(C=C2C(=C1)C(=O)OC3=C2C=CC(=C3)O)O"},
    
    # 14. Searsia pyroides
    {"Compound_Name": "Amentoflavone", "Plant_Source": "Searsia pyroides", "SMILES": "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)C4=C(C=CC(=C4)C5=CC(=O)C6=C(C=C(C=C6O5)O)O)O"},
    {"Compound_Name": "Fisetin", "Plant_Source": "Searsia pyroides", "SMILES": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(O2)C=C(C=C3)O)O)O)O"},
    
    # 15. Euclea crispa
    {"Compound_Name": "7-Methyljuglone", "Plant_Source": "Euclea crispa", "SMILES": "CC1=CC2=C(C(=C1)O)C(=O)C=CC2=O"},
    {"Compound_Name": "Euclein", "Plant_Source": "Euclea crispa", "SMILES": "CC1=CC2=C(C(=O)C=CC2=O)C(=C1)O"},
    
    # 16. Bridelia mollis
    {"Compound_Name": "Bridelionoside", "Plant_Source": "Bridelia mollis", "SMILES": "CC1C(C(C(C(O1)OCC2C(C(C(O2)O)O)O)O)O)O"},
    {"Compound_Name": "Isorhamnetin", "Plant_Source": "Bridelia mollis", "SMILES": "COC1=C(C=CC(=C1)C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O"},
    
    # 17. Vitex payos
    {"Compound_Name": "Vitexin", "Plant_Source": "Vitex payos", "SMILES": "C1=CC(=CC=C1C2=CC(=O)C3=C(C(=C(C=C3O2)O)C4C(C(C(C(O4)CO)O)O)O)O)O"},
    {"Compound_Name": "Casticin", "Plant_Source": "Vitex payos", "SMILES": "COC1=C(C=C(C=C1)C2=C(C(=O)C3=C(C(=C(C(=C3O2)OC)OC)O)OC)O)O"},
    
    # 18. Ficus ingens
    {"Compound_Name": "Rutin", "Plant_Source": "Ficus ingens", "SMILES": "CC1C(C(C(C(O1)OCC2C(C(C(C(O2)OC3=C(OC4=CC(=CC(=C4C3=O)O)O)C5=CC(=C(C=C5)O)O)O)O)O)O)O)O"},
    {"Compound_Name": "Lupeol", "Plant_Source": "Ficus ingens", "SMILES": "CC(=C)C1CCC2(C1C3CCC4C5(CCC(C(C5CCC4(C3(CC2)C)C)(C)C)O)C)C"},
    
    # 19. Psidium guajava
    {"Compound_Name": "Guaijaverin", "Plant_Source": "Psidium guajava", "SMILES": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)OC4C(C(C(O4)O)O)O)O)O"},
    {"Compound_Name": "Caryophyllene", "Plant_Source": "Psidium guajava", "SMILES": "CC1=CCCC(=C)C2CC1C2(C)C"},
    
    # 20. Eucalyptus camaldulensis
    {"Compound_Name": "Eucalyptol", "Plant_Source": "Eucalyptus camaldulensis", "SMILES": "CC1(C2CCC(O2)(CC1)C)C"},
    {"Compound_Name": "Macrocarpal A", "Plant_Source": "Eucalyptus camaldulensis", "SMILES": "CC1=C(C(=O)C2=C(C1=O)C(C3=C(O2)C(=O)C(=C(C3=O)C)C(C)C)C)C(C)C"},
    
    # 21. Citrus limon
    {"Compound_Name": "Limonene", "Plant_Source": "Citrus limon", "SMILES": "CC1=CCC(CC1)C(=C)C"},
    {"Compound_Name": "Hesperidin", "Plant_Source": "Citrus limon", "SMILES": "CC1C(C(C(C(O1)OCC2C(C(C(C(O2)OC3CC(=O)C4=C(C=C(C=C4O3)O)O)C5=CC(=C(C=C5)OC)O)O)O)O)O)O"},
    
    # Control
    {"Compound_Name": "PT70 (Control Inhibitor)", "Plant_Source": "Synthetic Control", "SMILES": "CCCCCC1=CC=C(C=C1)CC2=CC=C(O2)C(=O)O"}
]

# ==============================================================================
# 2. OPTIONAL AUTOMATED PUBCHEM MINING FUNCTION
# ==============================================================================
def fetch_pubchem_compounds(plant_name):
    """Mines additional CIDs directly from PubChem API for a plant organism."""
    compounds = []
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/concept/name/{plant_name}/cids/JSON"
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            cids = response.json().get('IdentifierList', {}).get('CID', [])[:5] # fetch top 5 per plant
            for cid in cids:
                smiles_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,IUPACName/JSON"
                s_res = requests.get(smiles_url, timeout=5)
                if s_res.status_code == 200:
                    props = s_res.json()['PropertyTable']['Properties'][0]
                    compounds.append({
                        "Compound_Name": props.get('IUPACName', f"CID_{cid}"),
                        "Plant_Source": plant_name,
                        "SMILES": props.get('CanonicalSMILES')
                    })
                time.sleep(0.2)
    except Exception:
        pass
    return compounds

# ==============================================================================
# 3. PROCESSING & VALIDATION
# ==============================================================================
print("=== Generating Zimbabwean Phytochemical Library ===")

# Build initial dataframe
df = pd.DataFrame(SEED_DATA)

# Validate SMILES using RDKit
valid_mask = []
for smi in df["SMILES"]:
    mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    valid_mask.append(mol is not None)

df = df[valid_mask].drop_duplicates(subset=["SMILES"])

# Output Directory Setup
output_dir = os.path.join("inha_virtual_screening", "data")
os.makedirs(output_dir, exist_ok=True)

csv_path = os.path.join(output_dir, "zimbabwe_phytochemicals.csv")
smi_path = os.path.join(output_dir, "zimbabwe_phytochemicals.smi")

# Save outputs
df.to_csv(csv_path, index=False)
df[["SMILES", "Compound_Name"]].to_csv(smi_path, sep="\t", index=False, header=False)

print(f"[+] Total Validated Compounds: {len(df)}")
print(f"[+] Total Unique Plant Species: {df[df['Plant_Source'] != 'Synthetic Control']['Plant_Source'].nunique()} / 21")
print(f"[+] Output CSV saved to : {csv_path}")
print(f"[+] Output SMI saved to : {smi_path}")