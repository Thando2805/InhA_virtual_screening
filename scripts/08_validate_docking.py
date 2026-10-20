import csv
import subprocess
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign


PROJECT_DIR = Path(__file__).resolve().parents[1]
STRUCTURE_FILE = PROJECT_DIR / "data" / "structures" / "2X22.pdb"
REFERENCE_SDF = PROJECT_DIR / "data" / "structures" / "TCU_ideal.sdf"
VINA_EXE = PROJECT_DIR / "bin" / "vina.exe"
OUTPUT_DIR = PROJECT_DIR / "validation" / "2X22_TCU_redocking"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def extract_crystal_ligand(pdb_lines: list[str]) -> Chem.Mol:
    ligand_lines = [
        line for line in pdb_lines
        if line.startswith("HETATM")
        and line[17:20] == "TCU"
        and line[21] == "A"
        and line[22:26].strip() == "1271"
    ]
    if not ligand_lines:
        raise ValueError("Could not find co-crystal TCU at chain A residue 1271 in 2X22.")

    ligand_serials = {int(line[6:11]) for line in ligand_lines}
    connectivity = []
    for line in pdb_lines:
        if not line.startswith("CONECT"):
            continue
        serials = [int(line[index:index + 5]) for index in range(6, len(line), 5) if line[index:index + 5].strip()]
        if serials and serials[0] in ligand_serials:
            connectivity.append(line)

    block = "\n".join(ligand_lines + connectivity + ["END"]) + "\n"
    ligand = Chem.MolFromPDBBlock(
        block,
        sanitize=False,
        removeHs=False,
        proximityBonding=False,
    )
    if ligand is None:
        raise ValueError("RDKit could not parse the crystal TCU coordinates and connectivity.")

    template = Chem.SDMolSupplier(str(REFERENCE_SDF), removeHs=False)[0]
    if template is None:
        raise ValueError(f"Could not read TCU reference chemistry from {REFERENCE_SDF}.")
    template = Chem.RemoveHs(template)
    ligand = AllChem.AssignBondOrdersFromTemplate(template, ligand)
    Chem.SanitizeMol(ligand)
    if ligand.GetNumAtoms() != 21:
        raise ValueError(f"Expected 21 heavy atoms in TCU; parsed {ligand.GetNumAtoms()}.")
    return ligand


def write_chain_a_receptor(pdb_lines: list[str], destination: Path) -> None:
    selected_atoms: dict[tuple[str, ...], tuple[float, str]] = {}
    for line in pdb_lines:
        if line.startswith("ATOM  ") and line[21] == "A":
            key = (
                line[21],
                line[22:27],
                line[12:16],
                line[76:78],
            )
            occupancy = float(line[54:60].strip() or 0.0)
            current = selected_atoms.get(key)
            altloc = line[16]
            if current is None or occupancy > current[0] or (
                occupancy == current[0] and altloc < current[1][16]
            ):
                selected_atoms[key] = (occupancy, line)

    receptor_lines = []
    for _, line in selected_atoms.values():
        if line[16] != " ":
            line = line[:16] + " " + line[17:]
        receptor_lines.append(line)
    receptor_lines.sort(key=lambda line: (int(line[22:26]), line[12:16]))
    if not receptor_lines:
        raise ValueError("No chain A protein atoms found in 2X22.")
    destination.write_text("\n".join(receptor_lines + ["TER", "END"]) + "\n", encoding="ascii")


def prepare_receptor(receptor_pdb: Path, receptor_pdbqt: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            "-m",
            "meeko.cli.mk_prepare_receptor",
            "--read_pdb",
            str(receptor_pdb),
            "--delete_bad_res",
            "--forgive_extra_bonds",
            "--write_pdbqt",
            str(receptor_pdbqt),
        ],
        check=True,
        cwd=PROJECT_DIR,
    )


def prepare_ligand(ligand: Chem.Mol, ligand_sdf: Path, ligand_pdbqt: Path) -> None:
    ligand_with_hydrogens = Chem.AddHs(ligand, addCoords=True)
    Chem.MolToMolFile(ligand_with_hydrogens, str(ligand_sdf))
    subprocess.run(
        [
            sys.executable,
            "-m",
            "meeko.cli.mk_prepare_ligand",
            "-i",
            str(ligand_sdf),
            "-o",
            str(ligand_pdbqt),
        ],
        check=True,
        cwd=PROJECT_DIR,
    )


def ligand_box(ligand: Chem.Mol) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    conformer = ligand.GetConformer()
    coordinates = [conformer.GetAtomPosition(index) for index in range(ligand.GetNumAtoms())]
    center = tuple(
        sum(point[axis] for point in coordinates) / len(coordinates)
        for axis in range(3)
    )
    size = tuple(
        max(20.0, max(point[axis] for point in coordinates) - min(point[axis] for point in coordinates) + 12.0)
        for axis in range(3)
    )
    return center, size


def parse_vina_heavy_atoms(pose_file: Path) -> Chem.Mol:
    pdb_records = []
    in_first_model = False
    serial = 1
    for line in pose_file.read_text(encoding="ascii").splitlines():
        if line.startswith("MODEL"):
            if in_first_model:
                break
            in_first_model = True
            continue
        if line.startswith("ENDMDL") and in_first_model:
            break
        if not in_first_model or not line.startswith(("ATOM", "HETATM")):
            continue

        fields = line.split()
        atom_type = fields[-1].strip()
        if atom_type.startswith("H"):
            continue
        element = "C" if atom_type == "A" else atom_type[0]
        x, y, z = (float(value) for value in fields[5:8])
        atom_name = line[12:16].strip() or f"{element}{serial}"
        pdb_records.append(
            f"HETATM{serial:5d} {atom_name:>4s} TCU A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}"
        )
        serial += 1

    docked = Chem.MolFromPDBBlock(
        "\n".join(pdb_records + ["END"]) + "\n",
        sanitize=False,
        removeHs=True,
        proximityBonding=True,
    )
    if docked is None:
        raise ValueError("Could not reconstruct heavy-atom connectivity from the top Vina pose.")
    template = Chem.RemoveHs(Chem.SDMolSupplier(str(REFERENCE_SDF), removeHs=False)[0])
    docked = AllChem.AssignBondOrdersFromTemplate(template, docked)
    Chem.SanitizeMol(docked)
    if docked.GetNumAtoms() != 21:
        raise ValueError(f"Expected 21 heavy atoms in top Vina pose; parsed {docked.GetNumAtoms()}.")
    return docked


def main() -> None:
    for required in (STRUCTURE_FILE, REFERENCE_SDF, VINA_EXE):
        if not required.is_file():
            raise FileNotFoundError(f"Required validation input not found: {required}")

    pdb_lines = STRUCTURE_FILE.read_text(encoding="ascii").splitlines()
    crystal_ligand = extract_crystal_ligand(pdb_lines)
    receptor_pdb = OUTPUT_DIR / "2X22_chainA_protein_only.pdb"
    receptor_pdbqt = OUTPUT_DIR / "2X22_chainA_protein_only.pdbqt"
    ligand_sdf = OUTPUT_DIR / "TCU_crystal_pose.sdf"
    ligand_pdbqt = OUTPUT_DIR / "TCU_crystal_pose.pdbqt"
    pose_file = OUTPUT_DIR / "TCU_redocked.pdbqt"
    result_file = OUTPUT_DIR / "redocking_validation.csv"

    write_chain_a_receptor(pdb_lines, receptor_pdb)
    prepare_receptor(receptor_pdb, receptor_pdbqt)
    prepare_ligand(crystal_ligand, ligand_sdf, ligand_pdbqt)

    center, size = ligand_box(crystal_ligand)
    vina_command = [
        str(VINA_EXE),
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}",
        "--size_y", f"{size[1]:.3f}",
        "--size_z", f"{size[2]:.3f}",
        "--exhaustiveness", "12",
        "--num_modes", "9",
        "--seed", "42",
        "--out", str(pose_file),
    ]
    if not pose_file.is_file():
        vina_run = subprocess.run(vina_command, check=True, capture_output=True, text=True)
        (OUTPUT_DIR / "vina.log").write_text(
            vina_run.stdout + "\n" + vina_run.stderr,
            encoding="utf-8",
        )
    else:
        print(f"Reusing existing docked pose: {pose_file}")

    affinity = None
    for line in pose_file.read_text(encoding="ascii").splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            affinity = float(line.split()[3])
            break
    if affinity is None:
        raise ValueError("Vina output did not include a binding affinity for MODEL 1.")

    docked_ligand = parse_vina_heavy_atoms(pose_file)
    rmsd = rdMolAlign.CalcRMS(crystal_ligand, docked_ligand)
    row = {
        "PDB": "2X22",
        "Crystal_Ligand": "TCU (PT70)",
        "Crystal_Ligand_Residue": "A:1271",
        "Receptor": "2X22 chain A protein only; TCU and NAD excluded to match project protocol",
        "Heavy_Atoms": crystal_ligand.GetNumAtoms(),
        "Vina_Affinity_kcal_mol": affinity,
        "Top_Pose_Symmetry_Aware_Heavy_Atom_RMSD_A": rmsd,
        "RMSD_Pass_At_2A": rmsd <= 2.0,
        "Box_Center_XYZ_A": ";".join(f"{value:.3f}" for value in center),
        "Box_Size_XYZ_A": ";".join(f"{value:.3f}" for value in size),
        "Exhaustiveness": 12,
        "Seed": 42,
    }
    with result_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)

    print(f"Top-pose affinity: {affinity:.3f} kcal/mol")
    print(f"Symmetry-aware heavy-atom RMSD: {rmsd:.3f} A")
    print(f"2.0 A redocking criterion: {'PASS' if rmsd <= 2.0 else 'FAIL'}")
    print(f"Results saved to {result_file}")
    print(f"Docked pose saved to {pose_file}")


if __name__ == "__main__":
    main()