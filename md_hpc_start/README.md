# HPC Molecular Dynamics Starter: InhA 2AQ8 Complexes

This folder is a practical bridge from the completed docking work to an HPC molecular dynamics run. It is not a finished MD system yet. The next required step is to build validated protein/cofactor/ligand topologies, then run minimization, equilibration, production, and analysis.

## What you already have

- Protein/crystal reference: `../galaxy_handoff_2AQ8/structures/2AQ8_full_crystal.pdb`
- Docking receptor/protein files: `../galaxy_handoff_2AQ8/structures/`
- Docked ligand poses: `../galaxy_handoff_2AQ8/ligands/docked_poses/`
- Docked SDF poses: `../galaxy_handoff_2AQ8/ligands/docked_poses_sdf/`
- Ligand SMILES: `../galaxy_handoff_2AQ8/results/original_2AQ8_ligand_smiles.csv`
- Docking rank summary: `../galaxy_handoff_2AQ8/results/final_docking_summary.csv`

## Recommended pilot systems

Start with only two systems before scaling up:

1. `Lantadene_A` as the best docking hit.
2. `PT70_Control_Inhibitor` as the control.

After both are stable through minimization, NVT, NPT, and a short production test, expand to Rutin, Guaijaverin, Quercitrin, and 7-Methyljuglone.

## Critical choice before running MD

For InhA, decide whether to keep NAD from the crystal structure. The full crystal file contains NAD as residue `NAI` in chain A. The docking receptor removed it, but that docking decision should not automatically carry into MD.

Recommended starting hypothesis: keep NAD/cofactor for the biological enzyme system, unless your supervisor specifically wants apo InhA.

## Suggested GROMACS preparation route

Use this if your HPC has GROMACS.

1. Prepare protein from `2AQ8_full_crystal.pdb`.
2. Keep or remove NAD intentionally.
3. Generate ligand parameters from the verified ligand structure/SMILES, not from PDBQT as if it were a topology.
4. Place the ligand using the docked SDF/PDBQT pose after checking atom mapping and protonation.
5. Combine protein, optional NAD, ligand, water, and ions.
6. Run:
   - energy minimization
   - NVT equilibration
   - NPT equilibration
   - short production test
   - full production replicas

## Files in this starter pack

- `mdp/minim.mdp`: energy minimization settings
- `mdp/nvt.mdp`: constant volume equilibration
- `mdp/npt.mdp`: constant pressure equilibration
- `mdp/md_10ns_test.mdp`: short test production run
- `slurm/gromacs_pipeline_template.sbatch`: HPC job template

## Before submitting to HPC

Confirm these details on the cluster:

- Scheduler: SLURM, PBS, or another system
- GROMACS module name, for example `gromacs/2023`
- CPU/GPU availability
- Walltime and partition names
- Whether ligand topology generation tools are installed, such as ACPYPE, CGenFF, OpenFF, AmberTools, or GAFF workflow tools

## Minimum analysis after each pilot

Check these before trusting the simulation:

- Potential energy after minimization
- Temperature stability during NVT
- Pressure/density stability during NPT
- Protein backbone RMSD
- Ligand RMSD in the binding pocket
- Protein-ligand hydrogen bonds/contacts
- Whether ligand leaves the pocket

