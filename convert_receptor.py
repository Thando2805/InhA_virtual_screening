from pathlib import Path
import subprocess


project_dir = Path(__file__).resolve().parent
pdb_file = project_dir / "data" / "structures" / "2aq8_chainA_clean.pdb"
output_base = project_dir / "2aq8_receptor"

subprocess.run(
    [
        "mk_prepare_receptor.exe",
        "--read_pdb",
        str(pdb_file),
        "--delete_residues",
        "A:270,447",
        "--delete_bad_res",
        "--forgive_extra_bonds",
        "-o",
        str(output_base),
        "-p",
    ],
    check=True,
)