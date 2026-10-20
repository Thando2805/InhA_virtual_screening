# HPC Build Guide

This package contains docking poses and MD parameter templates. Before running GROMACS, each selected ligand-protein system must be converted into a complete, force-field-compatible system.

## Selected Systems

Build in this order:

1. `systems/Lantadene_A`
2. `systems/PT70_Control_Inhibitor`
3. `systems/Rutin`

The SLURM script expects each system directory to contain:

```text
complex.gro
topol.top
index.ndx
```

## Shared Inputs

Protein:

```text
protein/2aq8_chainA_clean.pdb
```

Docked ligand poses:

```text
ligands/Lantadene_A/Lantadene_A_docked.sdf
ligands/PT70_Control_Inhibitor/PT70_Control_Inhibitor_docked.sdf
ligands/Rutin/Rutin_docked.sdf
```

SMILES are stored in:

```text
selected_compounds.csv
```

## Recommended Route A: CHARMM36m + CGenFF

Use this if the HPC has CHARMM36m and CGenFF support.

General workflow:

```bash
cd systems/Lantadene_A
cp ../../protein/2aq8_chainA_clean.pdb protein.pdb
cp ../../ligands/Lantadene_A/Lantadene_A_docked.sdf ligand_docked.sdf

# Generate protein topology with CHARMM36m.
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -p topol.top -i posre.itp -water tip3p

# Generate ligand parameters from verified ligand chemistry using CGenFF.
# Follow the local CGenFF/CHARMM-GUI workflow available on your HPC.

# Merge protein and ligand coordinates into complex.gro.
# Include ligand .itp files in topol.top.
# Solvate, add ions, and create index.ndx.
```

Use the same force-field route for all three systems.

## Recommended Route B: AMBER Protein + GAFF/ACPYPE Ligands

Use this if the HPC has AmberTools, ACPYPE, or a lab-standard AMBER-to-GROMACS workflow.

General workflow:

```bash
cd systems/Lantadene_A
cp ../../protein/2aq8_chainA_clean.pdb protein.pdb
cp ../../ligands/Lantadene_A/Lantadene_A_docked.sdf ligand_docked.sdf

# Generate protein topology with a compatible AMBER force field.
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -p topol.top -i posre.itp -water tip3p

# Generate ligand parameters from SDF/SMILES using GAFF/ACPYPE.
# Confirm ligand charge and protonation before accepting the topology.

# Merge coordinates, include ligand topology, solvate, add ions, and create index.ndx.
```

## Index Group Requirement

The included MDP files use:

```text
tc-grps = Protein_LIG Water_and_ions
```

Create an index file with those exact groups, or edit the MDP files to match your generated group names.

Typical `make_ndx` idea:

```bash
gmx make_ndx -f complex.gro -o index.ndx
```

Create:

```text
Protein_LIG
Water_and_ions
```

## Run Template

After `complex.gro`, `topol.top`, and `index.ndx` exist:

```bash
cd systems/Lantadene_A
cp ../../slurm/gromacs_pipeline_template.sbatch .
sbatch gromacs_pipeline_template.sbatch
```

Repeat for the PT70 control only after the Lantadene_A pilot setup runs cleanly.

## Pre-Submission Checklist

- Confirm ligand protonation state and formal charge.
- Confirm ligand topology atom names match ligand coordinates.
- Confirm ligand pose is still aligned with the docked binding site.
- Confirm `topol.top` includes all ligand `.itp` files.
- Confirm the system is solvated and neutralized with about 0.15 M NaCl.
- Confirm `index.ndx` contains `Protein_LIG` and `Water_and_ions`.
- Run `gmx grompp` interactively once before submitting a long job.

