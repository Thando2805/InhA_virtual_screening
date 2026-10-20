import os
import pandas as pd
from chembl_webresource_client.new_client import new_client

os.makedirs("data/raw", exist_ok=True)

print("Connecting to ChEMBL database...")
target_api = new_client.target
targets = target_api.filter(target_components__accession='P9WGR1')
target_chembl_id = targets[0]['target_chembl_id']
print(f"Target ChEMBL ID: {target_chembl_id}")

activity_api = new_client.activity

# Keep the source dataset consistent with the IC50-based preprocessing.
res_ic50 = activity_api.filter(target_chembl_id=target_chembl_id, standard_type="IC50")
df = pd.DataFrame(list(res_ic50))
df = df.dropna(subset=['standard_value', 'canonical_smiles'])

output_path = "data/raw/inhA_chembl_ic50_raw.csv"
df.to_csv(output_path, index=False)
print(f"Downloaded {len(df)} total bioactivity records (IC50 & Ki) to {output_path}")