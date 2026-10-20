import csv
import subprocess
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


PROJECT_DIR = Path(__file__).resolve().parents[1]
HITS_FILE = PROJECT_DIR / "data" / "processed" / "zimbabwean_phytochemical_hits.csv"
RECEPTOR = PROJECT_DIR / "validation" / "2X22_TCU_redocking" / "2X22_chainA_protein_only.pdbqt"
TCU_REFERENCE = PROJECT_DIR / "data" / "structures" / "TCU_ideal.sdf"
TCU_CRYSTAL_PDBQT = PROJECT_DIR / "validation" / "2X22_TCU_redocking" / "TCU_crystal_pose.pdbqt"
VINA = PROJECT_DIR / "bin" / "vina.exe"
OUTPUT_DIR = PROJECT_DIR / "validation" / "exploratory_2X22"
LIGAND_DIR = OUTPUT_DIR / "ligands"
POSE_DIR = OUTPUT_DIR / "poses"
SUMMARY = OUTPUT_DIR / "exploratory_docking_summary.csv"

# Use one fixed search box centered on the validated TCU binding site.
BOX_CENTER = (-16.190, -35.898, 15.695)
BOX_SIZE = (21.484, 20.000, 20.000)
EXHAUSTIVENESS = 12
NUM_MODES = 9
SEED = 42


def prepare_candidate(name: str, smiles: str, output_pdbqt: Path, sdf_path: Path) -> int:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES for {name}")

    molecule = Chem.AddHs(molecule)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(molecule, params) != 0:
        raise RuntimeError(f"Could not generate a 3D conformer for {name}")
    AllChem.MMFFOptimizeMolecule(molecule)
    heavy_atoms = molecule.GetNumHeavyAtoms()
    Chem.MolToMolFile(molecule, str(sdf_path))
    subprocess.run(
        [
            sys.executable,
            "-m",
            "meeko.cli.mk_prepare_ligand",
            "-i",
            str(sdf_path),
            "-o",
            str(output_pdbqt),
        ],
        check=True,
        cwd=PROJECT_DIR,
    )
    return heavy_atoms


def run_docking(name: str, ligand_file: Path, output_file: Path) -> tuple[float, Path]:
    command = [
        str(VINA),
        "--receptor", str(RECEPTOR),
        "--ligand", str(ligand_file),
        "--center_x", str(BOX_CENTER[0]),
        "--center_y", str(BOX_CENTER[1]),
        "--center_z", str(BOX_CENTER[2]),
        "--size_x", str(BOX_SIZE[0]),
        "--size_y", str(BOX_SIZE[1]),
        "--size_z", str(BOX_SIZE[2]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--seed", str(SEED),
        "--out", str(output_file),
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    log_path = OUTPUT_DIR / f"{name}_vina.log"
    log_path.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")

    affinity = None
    for line in output_file.read_text(encoding="ascii").splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            affinity = float(line.split()[3])
            break
    if affinity is None:
        raise RuntimeError(f"No Vina affinity was found in {output_file}")
    return affinity, log_path


def main() -> None:
    for path in (HITS_FILE, RECEPTOR, TCU_REFERENCE, TCU_CRYSTAL_PDBQT, VINA):
        if not path.is_file():
            raise FileNotFoundError(f"Required input is missing: {path}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LIGAND_DIR.mkdir(parents=True, exist_ok=True)
    POSE_DIR.mkdir(parents=True, exist_ok=True)

    hits = pd.read_csv(HITS_FILE)
    natural_candidates = hits[hits["Plant_Source"] != "Synthetic Control"].head(5)
    if len(natural_candidates) != 5:
        raise ValueError(f"Expected five ranked natural candidates; found {len(natural_candidates)}")

    rows = []
    tcu = Chem.SDMolSupplier(str(TCU_REFERENCE), removeHs=False)[0]
    if tcu is None:
        raise ValueError(f"Could not read TCU structure from {TCU_REFERENCE}")
    tcu_heavy_atoms = Chem.RemoveHs(tcu).GetNumAtoms()
    control_pose = POSE_DIR / "TCU_crystal_control_docked.pdbqt"
    control_affinity, _ = run_docking("TCU_crystal_control", TCU_CRYSTAL_PDBQT, control_pose)
    rows.append({
        "Compound": "TCU (co-crystal PT70 control)",
        "Plant_Source_or_Role": "Positive control from 2X22",
        "Current_ML_Score": "",
        "Model_Threshold_Result": "Control; not scored in candidate panel",
        "Vina_Binding_Affinity_kcal_mol": control_affinity,
        "Heavy_Atom_Count": tcu_heavy_atoms,
        "Docking_Ligand_Efficiency_kcal_mol_per_heavy_atom": -control_affinity / tcu_heavy_atoms,
        "Docking_Set": "Exploratory; 2X22 chain A protein-only receptor",
        "Pose_File": str(control_pose.relative_to(PROJECT_DIR)),
    })
    print(f"Docked TCU crystal control: {control_affinity:.3f} kcal/mol")

    for candidate in natural_candidates.itertuples(index=False):
        safe_name = "".join(char if char.isalnum() or char in "-_." else "_" for char in candidate.Compound_Name)
        ligand_pdbqt = LIGAND_DIR / f"{safe_name}.pdbqt"
        ligand_sdf = LIGAND_DIR / f"{safe_name}.sdf"
        pose_file = POSE_DIR / f"{safe_name}_docked.pdbqt"
        heavy_atoms = prepare_candidate(candidate.Compound_Name, candidate.SMILES, ligand_pdbqt, ligand_sdf)
        affinity, _ = run_docking(safe_name, ligand_pdbqt, pose_file)
        rows.append({
            "Compound": candidate.Compound_Name,
            "Plant_Source_or_Role": candidate.Plant_Source,
            "Current_ML_Score": float(candidate.Inhibition_Probability),
            "Model_Threshold_Result": "Below 0.50; exploratory only",
            "Vina_Binding_Affinity_kcal_mol": affinity,
            "Heavy_Atom_Count": heavy_atoms,
            "Docking_Ligand_Efficiency_kcal_mol_per_heavy_atom": -affinity / heavy_atoms,
            "Docking_Set": "Exploratory; 2X22 chain A protein-only receptor",
            "Pose_File": str(pose_file.relative_to(PROJECT_DIR)),
        })
        print(f"Docked {candidate.Compound_Name}: {affinity:.3f} kcal/mol; ML score {candidate.Inhibition_Probability:.4f}")

    results = pd.DataFrame(rows).sort_values(
        "Vina_Binding_Affinity_kcal_mol", ascending=True
    )
    results["Docking_Box_Center_XYZ_A"] = ";".join(f"{value:.3f}" for value in BOX_CENTER)
    results["Docking_Box_Size_XYZ_A"] = ";".join(f"{value:.3f}" for value in BOX_SIZE)
    results["Exhaustiveness"] = EXHAUSTIVENESS
    results["Num_Modes"] = NUM_MODES
    results["Random_Seed"] = SEED
    results.to_csv(SUMMARY, index=False)
    print(f"Exploratory docking results saved to {SUMMARY}")


if __name__ == "__main__":
    main()
