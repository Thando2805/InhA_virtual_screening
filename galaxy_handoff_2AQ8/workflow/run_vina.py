import csv
import re
import subprocess
from pathlib import Path


project_dir = Path(__file__).resolve().parent
vina_exe = project_dir / "bin" / "vina.exe"
receptor_file = project_dir / "2aq8_receptor.pdbqt"
input_dir = project_dir / "ml_selected_pdbqt"
output_dir = project_dir / "docking_results"
summary_file = project_dir / "final_docking_summary.csv"
box_center = (2.2, -32.2, 14.0)
box_size = (22.0, 22.0, 22.0)

for required_file in (vina_exe, receptor_file):
    if not required_file.is_file():
        raise FileNotFoundError(f"Required file not found: {required_file}")

ligand_files = sorted(input_dir.glob("*.pdbqt"))
if not ligand_files:
    raise FileNotFoundError(f"No ligand PDBQT files found in {input_dir}")

output_dir.mkdir(parents=True, exist_ok=True)
results = []

print(f"--- Docking {len(ligand_files)} compounds against InhA (2AQ8) ---\n")

for ligand_file in ligand_files:
    ligand_name = ligand_file.stem
    pose_file = output_dir / f"{ligand_name}_docked.pdbqt"
    command = [
        str(vina_exe),
        "--receptor", str(receptor_file),
        "--ligand", str(ligand_file),
        "--center_x", str(box_center[0]),
        "--center_y", str(box_center[1]),
        "--center_z", str(box_center[2]),
        "--size_x", str(box_size[0]),
        "--size_y", str(box_size[1]),
        "--size_z", str(box_size[2]),
        "--exhaustiveness", "12",
        "--num_modes", "5",
        "--out", str(pose_file),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)

    pose_text = pose_file.read_text(encoding="ascii")
    match = re.search(r"^REMARK VINA RESULT:\s+([-+]?\d+(?:\.\d+)?)", pose_text, re.MULTILINE)
    if match is None:
        raise RuntimeError(f"Could not read Vina affinity from {pose_file}")

    affinity = float(match.group(1))
    results.append((ligand_name, affinity))
    print(f"Docked {ligand_name}: {affinity:.2f} kcal/mol")

results.sort(key=lambda result: result[1])
with summary_file.open("w", newline="", encoding="utf-8") as summary:
    writer = csv.writer(summary)
    writer.writerow(("Compound", "Vina_Binding_Affinity_kcal_mol"))
    writer.writerows(results)

print(f"\nResults saved to {summary_file}")