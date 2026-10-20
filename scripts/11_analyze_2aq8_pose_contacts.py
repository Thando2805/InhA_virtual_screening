import csv
from collections import defaultdict
from pathlib import Path

import numpy as np


PROJECT_DIR = Path(__file__).resolve().parents[1]
RECEPTOR_FILE = PROJECT_DIR / "2aq8_receptor.pdbqt"
SUMMARY_FILE = PROJECT_DIR / "final_docking_summary.csv"
POSE_DIR = PROJECT_DIR / "docking_results"
OUTPUT_DIR = PROJECT_DIR / "validation" / "original_2AQ8_contacts"
RESIDUE_SUMMARY = OUTPUT_DIR / "2aq8_pose_contact_summary.csv"
PAIR_DETAILS = OUTPUT_DIR / "2aq8_close_contact_pairs.csv"

CONTACT_CUTOFF_A = 4.0
POLAR_CUTOFF_A = 3.5


def parse_pdbqt_atoms(path: Path, ligand: bool = False) -> list[dict[str, object]]:
    atoms = []
    in_first_model = not ligand
    for line in path.read_text(encoding="ascii").splitlines():
        if ligand and line.startswith("MODEL"):
            in_first_model = line.split()[1] == "1"
            continue
        if ligand and line.startswith("ENDMDL") and in_first_model:
            break
        if not in_first_model or not line.startswith(("ATOM", "HETATM")):
            continue

        fields = line.split()
        atom_type = fields[-1].strip()
        if atom_type.startswith("H"):
            continue
        residue_name = fields[3]
        if ligand:
            chain = ""
            residue_number = fields[4]
            coordinate_fields = fields[5:8]
        else:
            chain = fields[4]
            residue_number = fields[5]
            coordinate_fields = fields[6:9]
        coordinates = np.asarray([float(value) for value in coordinate_fields], dtype=float)
        atoms.append({
            "atom_name": fields[2],
            "residue_name": residue_name,
            "chain": chain,
            "residue_number": residue_number,
            "atom_type": atom_type,
            "element": "C" if atom_type == "A" else atom_type[0],
            "coordinates": coordinates,
        })
    return atoms


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    receptor_atoms = parse_pdbqt_atoms(RECEPTOR_FILE)
    if not receptor_atoms:
        raise ValueError(f"No heavy atoms parsed from receptor: {RECEPTOR_FILE}")

    with SUMMARY_FILE.open(newline="", encoding="utf-8") as handle:
        compounds = list(csv.DictReader(handle))

    summary_rows = []
    pair_rows = []
    receptor_coords = np.vstack([atom["coordinates"] for atom in receptor_atoms])
    receptor_polar = np.asarray([atom["element"] in {"N", "O", "S"} for atom in receptor_atoms])

    for compound in compounds:
        name = compound["Compound"]
        pose_path = POSE_DIR / f"{name}_docked.pdbqt"
        ligand_atoms = parse_pdbqt_atoms(pose_path, ligand=True)
        if not ligand_atoms:
            raise ValueError(f"No top-pose heavy atoms parsed for {name} from {pose_path}")

        ligand_coords = np.vstack([atom["coordinates"] for atom in ligand_atoms])
        distances = np.linalg.norm(ligand_coords[:, None, :] - receptor_coords[None, :, :], axis=2)
        close_indices = np.argwhere(distances <= CONTACT_CUTOFF_A)
        polar_indices = np.argwhere(
            (distances <= POLAR_CUTOFF_A)
            & np.asarray([atom["element"] in {"N", "O", "S"} for atom in ligand_atoms])[:, None]
            & receptor_polar[None, :]
        )

        residues = defaultdict(lambda: {"min_distance": float("inf"), "atom_pairs": 0})
        for ligand_index, receptor_index in close_indices:
            receptor_atom = receptor_atoms[receptor_index]
            residue_key = (
                f"{receptor_atom['residue_name']}{receptor_atom['residue_number']}"
                f"{receptor_atom['chain']}"
            )
            residue = residues[residue_key]
            residue["min_distance"] = min(residue["min_distance"], float(distances[ligand_index, receptor_index]))
            residue["atom_pairs"] += 1
            pair_rows.append({
                "Compound": name,
                "Ligand_Atom": ligand_atoms[ligand_index]["atom_name"],
                "Ligand_Atom_Type": ligand_atoms[ligand_index]["atom_type"],
                "Protein_Residue": residue_key,
                "Protein_Atom": receptor_atom["atom_name"],
                "Protein_Atom_Type": receptor_atom["atom_type"],
                "Distance_A": float(distances[ligand_index, receptor_index]),
                "Close_Polar_Geometry_Candidate": (
                    distances[ligand_index, receptor_index] <= POLAR_CUTOFF_A
                    and ligand_atoms[ligand_index]["element"] in {"N", "O", "S"}
                    and receptor_atom["element"] in {"N", "O", "S"}
                ),
            })

        ordered_residues = sorted(residues.items(), key=lambda item: item[1]["min_distance"])
        summary_rows.append({
            "Compound": name,
            "Vina_Binding_Affinity_kcal_mol": compound["Vina_Binding_Affinity_kcal_mol"],
            "Heavy_Atoms": len(ligand_atoms),
            "Protein_Residues_with_heavy_atom_contact_within_4A": "; ".join(
                f"{residue} ({details['min_distance']:.2f} A)"
                for residue, details in ordered_residues
            ),
            "Residue_Count_within_4A": len(ordered_residues),
            "Heavy_Atom_Pairs_within_4A": len(close_indices),
            "Close_Polar_Geometry_Pairs_within_3_5A": len(polar_indices),
            "Interpretation": "Geometric contacts only; polar pairs are not hydrogen-bond assignments without angle/protonation review.",
        })

    with RESIDUE_SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    with PAIR_DETAILS.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)

    for row in summary_rows:
        print(f"{row['Compound']}: {row['Residue_Count_within_4A']} contacting residues; "
              f"{row['Heavy_Atom_Pairs_within_4A']} atom pairs <={CONTACT_CUTOFF_A:.1f} A; "
              f"{row['Close_Polar_Geometry_Pairs_within_3_5A']} close polar pairs")
    print(f"Residue summary: {RESIDUE_SUMMARY}")
    print(f"Atom-pair details: {PAIR_DETAILS}")


if __name__ == "__main__":
    main()
