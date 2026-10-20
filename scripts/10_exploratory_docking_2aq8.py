import subprocess
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


PROJECT_DIR = Path(__file__).resolve().parents[1]
HITS_FILE = PROJECT_DIR / "data" / "processed" / "zimbabwean_phytochemical_hits.csv"
RECEPTOR = PROJECT_DIR / "2aq8_receptor.pdbqt"
TCU_REFERENCE = PROJECT_DIR / "data" / "structures" / "TCU_ideal.sdf"
VINA = PROJECT_DIR / "bin" / "vina.exe"
OUTPUT_DIR = PROJECT_DIR / "validation" / "exploratory_2AQ8"
LIGAND_DIR = OUTPUT_DIR / "ligands"
POSE_DIR = OUTPUT_DIR / "poses"
SUMMARY = OUTPUT_DIR / "exploratory_docking_summary.csv"

# Match the search box and settings used by the original 2AQ8 run_vina.py.
BOX_CENTER = (2.2, -32.2, 14.0)
BOX_SIZE = (22.0, 22.0, 22.0)
EXHAUSTIVENESS = 12
NUM_MODES = 9
SEED = 42


def prepare_smiles_ligand(name: str, smiles: str, output_pdbqt: Path, output_sdf: Path) -> int:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES for {name}")
    molecule = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = SEED
    if AllChem.EmbedMolecule(molecule, parameters) != 0:
        raise RuntimeError(f"Could not create a 3D conformer for {name}")
    AllChem.MMFFOptimizeMolecule(molecule)
    heavy_atoms = molecule.GetNumHeavyAtoms()
    Chem.MolToMolFile(molecule, str(output_sdf))
    subprocess.run(
        [
            sys.executable,
            "-m",
            "meeko.cli.mk_prepare_ligand",
            "-i",
            str(output_sdf),
            "-o",
            str(output_pdbqt),
        ],
        check=True,
        cwd=PROJECT_DIR,
    )
    return heavy_atoms


def prepare_tcu_control(output_pdbqt: Path) -> int:
    reference = Chem.SDMolSupplier(str(TCU_REFERENCE), removeHs=False)[0]
    if reference is None:
        raise ValueError(f"Could not load TCU/PT70 reference from {TCU_REFERENCE}")
    heavy_atoms = Chem.RemoveHs(reference).GetNumAtoms()
    subprocess.run(
        [
            sys.executable,
            "-m",
            "meeko.cli.mk_prepare_ligand",
            "-i",
            str(TCU_REFERENCE),
            "-o",
            str(output_pdbqt),
        ],
        check=True,
        cwd=PROJECT_DIR,
    )
    return heavy_atoms


def dock(name: str, ligand_file: Path, pose_file: Path) -> float:
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
        "--out", str(pose_file),
    ]
    run = subprocess.run(command, check=True, capture_output=True, text=True)
    (OUTPUT_DIR / f"{name}_vina.log").write_text(run.stdout + "\n" + run.stderr, encoding="utf-8")
    for line in pose_file.read_text(encoding="ascii").splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    raise RuntimeError(f"No Vina score found in {pose_file}")


def make_result(name: str, role: str, ml_score: object, threshold: str, affinity: float,
                heavy_atoms: int, pose_file: Path) -> dict[str, object]:
    return {
        "Compound": name,
        "Plant_Source_or_Role": role,
        "Current_ML_Score": ml_score,
        "Model_Threshold_Result": threshold,
        "Vina_Binding_Affinity_kcal_mol": affinity,
        "Heavy_Atom_Count": heavy_atoms,
        "Docking_Ligand_Efficiency_kcal_mol_per_heavy_atom": -affinity / heavy_atoms,
        "Docking_Set": "Exploratory; existing 2AQ8 receptor; original run_vina.py box",
        "Pose_File": str(pose_file.relative_to(PROJECT_DIR)),
        "Docking_Box_Center_XYZ_A": ";".join(f"{value:.3f}" for value in BOX_CENTER),
        "Docking_Box_Size_XYZ_A": ";".join(f"{value:.3f}" for value in BOX_SIZE),
        "Exhaustiveness": EXHAUSTIVENESS,
        "Num_Modes": NUM_MODES,
        "Random_Seed": SEED,
    }


def main() -> None:
    for path in (HITS_FILE, RECEPTOR, TCU_REFERENCE, VINA):
        if not path.is_file():
            raise FileNotFoundError(f"Required input is missing: {path}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LIGAND_DIR.mkdir(parents=True, exist_ok=True)
    POSE_DIR.mkdir(parents=True, exist_ok=True)

    hits = pd.read_csv(HITS_FILE)
    natural = hits[hits["Plant_Source"] != "Synthetic Control"].head(5)
    if len(natural) != 5:
        raise ValueError(f"Expected five ranked natural candidates, found {len(natural)}")

    results = []
    tcu_pdbqt = LIGAND_DIR / "TCU_PT70_control.pdbqt"
    tcu_pose = POSE_DIR / "TCU_PT70_control_docked.pdbqt"
    tcu_atoms = prepare_tcu_control(tcu_pdbqt)
    tcu_affinity = dock("TCU_PT70_control", tcu_pdbqt, tcu_pose)
    results.append(make_result(
        "TCU (co-crystal PT70 control)",
        "Positive control from 2X22; docked to 2AQ8",
        "",
        "Control; not a natural candidate",
        tcu_affinity,
        tcu_atoms,
        tcu_pose,
    ))
    print(f"Docked TCU/PT70 control to 2AQ8: {tcu_affinity:.3f} kcal/mol")

    for candidate in natural.itertuples(index=False):
        safe_name = "".join(char if char.isalnum() or char in "-_." else "_" for char in candidate.Compound_Name)
        ligand_file = LIGAND_DIR / f"{safe_name}.pdbqt"
        ligand_sdf = LIGAND_DIR / f"{safe_name}.sdf"
        pose_file = POSE_DIR / f"{safe_name}_docked.pdbqt"
        atom_count = prepare_smiles_ligand(candidate.Compound_Name, candidate.SMILES, ligand_file, ligand_sdf)
        affinity = dock(safe_name, ligand_file, pose_file)
        results.append(make_result(
            candidate.Compound_Name,
            candidate.Plant_Source,
            float(candidate.Inhibition_Probability),
            "Below 0.50; exploratory only",
            affinity,
            atom_count,
            pose_file,
        ))
        print(f"Docked {candidate.Compound_Name} to 2AQ8: {affinity:.3f} kcal/mol; "
              f"ML score {candidate.Inhibition_Probability:.4f}")

    results_df = pd.DataFrame(results).sort_values("Vina_Binding_Affinity_kcal_mol")
    results_df.to_csv(SUMMARY, index=False)
    print(f"2AQ8 exploratory results saved to {SUMMARY}")


if __name__ == "__main__":
    main()
