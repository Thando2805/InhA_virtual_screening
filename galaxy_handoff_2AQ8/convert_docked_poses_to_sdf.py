import csv
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Geometry import Point3D


HANDOFF_DIR = Path(__file__).resolve().parent
POSE_DIR = HANDOFF_DIR / "ligands" / "docked_poses"
SMILES_FILE = HANDOFF_DIR / "results" / "original_2AQ8_ligand_smiles.csv"
OUTPUT_DIR = HANDOFF_DIR / "ligands" / "docked_poses_sdf"


def read_top_pose(path: Path) -> tuple[dict[int, int], dict[int, int], dict[int, tuple[float, float, float]], dict[int, str], float | None]:
    smiles_indices = {}
    hydrogen_parents = {}
    coordinates = {}
    atom_types = {}
    affinity = None
    in_first_model = False

    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith("MODEL"):
            if in_first_model:
                break
            in_first_model = line.split()[1] == "1"
            continue
        if line.startswith("ENDMDL") and in_first_model:
            break
        if not in_first_model:
            continue
        if line.startswith("REMARK SMILES IDX"):
            values = [int(value) for value in line.split()[3:]]
            smiles_indices.update(
                (values[index + 1], values[index] - 1)
                for index in range(0, len(values), 2)
            )
        elif line.startswith("REMARK H PARENT"):
            values = [int(value) for value in line.split()[3:]]
            hydrogen_parents.update(
                (values[index + 1], values[index] - 1)
                for index in range(0, len(values), 2)
            )
        elif line.startswith("REMARK VINA RESULT:"):
            affinity = float(line.split()[3])
        elif line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            coordinates[serial] = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
            atom_types[serial] = line.split()[-1]

    if not coordinates or not smiles_indices:
        raise ValueError(f"Missing top-pose atoms or atom mapping in {path}")
    return smiles_indices, hydrogen_parents, coordinates, atom_types, affinity


def pdbqt_atomic_number(atom_type: str) -> int:
    if atom_type == "A":
        return 6
    symbol = atom_type[:1].upper()
    atomic_number = Chem.GetPeriodicTable().GetAtomicNumber(symbol)
    if not atomic_number:
        raise ValueError(f"Unsupported PDBQT atom type {atom_type!r} in {path}")
    return atomic_number


def convert(compound: str, smiles: str) -> Path:
    pose_path = POSE_DIR / f"{compound}_docked.pdbqt"
    smiles_indices, hydrogen_parents, coordinates, atom_types, affinity = read_top_pose(pose_path)
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Could not parse SMILES for {compound}")
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    expected_cip = {
        atom.GetIdx(): atom.GetProp("_CIPCode")
        for atom in molecule.GetAtoms()
        if atom.HasProp("_CIPCode")
    }

    heavy_coordinates = {}
    for serial, smiles_index in smiles_indices.items():
        if serial not in coordinates:
            raise ValueError(f"Mapped atom {serial} is absent from {pose_path}")
        if molecule.GetAtomWithIdx(smiles_index).GetAtomicNum() != pdbqt_atomic_number(atom_types[serial]):
            raise ValueError(f"Element mismatch for atom {serial} in {pose_path}")
        heavy_coordinates[smiles_index] = coordinates[serial]
    if set(heavy_coordinates) != set(range(molecule.GetNumAtoms())):
        raise ValueError(f"Heavy-atom mapping is incomplete for {compound}")

    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for atom_index, xyz in heavy_coordinates.items():
        conformer.SetAtomPosition(atom_index, Point3D(*xyz))
    molecule.AddConformer(conformer)
    molecule = Chem.AddHs(molecule, addCoords=True)

    hydrogens_by_parent = defaultdict(list)
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 1:
            hydrogens_by_parent[atom.GetNeighbors()[0].GetIdx()].append(atom.GetIdx())
    pose_hydrogens_by_parent = defaultdict(list)
    for hydrogen_serial, parent_index in hydrogen_parents.items():
        if hydrogen_serial not in coordinates or parent_index not in heavy_coordinates:
            raise ValueError(f"Invalid hydrogen-parent mapping in {pose_path}")
        pose_hydrogens_by_parent[parent_index].append(hydrogen_serial)
    for parent_index, hydrogen_indices in hydrogens_by_parent.items():
        pose_serials = sorted(pose_hydrogens_by_parent[parent_index])
        if len(pose_serials) > len(hydrogen_indices):
            raise ValueError(f"Hydrogen count mismatch for {compound}, atom {parent_index + 1}")
        for atom_index, serial in zip(hydrogen_indices, pose_serials):
            molecule.GetConformer().SetAtomPosition(atom_index, Point3D(*coordinates[serial]))

    pose_stereo = Chem.Mol(molecule)
    Chem.AssignStereochemistryFrom3D(pose_stereo, replaceExistingTags=True)
    Chem.AssignStereochemistry(pose_stereo, cleanIt=True, force=True)
    observed_cip = {
        atom.GetIdx(): atom.GetProp("_CIPCode")
        for atom in pose_stereo.GetAtoms()
        if atom.HasProp("_CIPCode")
    }
    if any(observed_cip.get(index) != cip for index, cip in expected_cip.items()):
        raise ValueError(f"Pose stereochemistry differs from source SMILES for {compound}: {expected_cip} != {observed_cip}")

    molecule.SetProp("_Name", compound)
    molecule.SetProp("Pose_Model", "1")
    if affinity is not None:
        molecule.SetProp("Vina_Binding_Affinity_kcal_mol", f"{affinity:.3f}")
    Chem.SanitizeMol(molecule)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / f"{compound}_docked.sdf"
    Chem.MolToMolFile(molecule, str(output_path))
    check = Chem.SDMolSupplier(str(output_path), removeHs=False)[0]
    if check is None or check.GetNumAtoms() != molecule.GetNumAtoms():
        raise ValueError(f"SDF validation failed for {output_path}")
    return output_path


def main() -> None:
    with SMILES_FILE.open(newline="", encoding="utf-8-sig") as handle:
        compounds = list(csv.DictReader(handle))
    for row in compounds:
        output_path = convert(row["Compound"], row["SMILES"])
        print(output_path)


if __name__ == "__main__":
    main()