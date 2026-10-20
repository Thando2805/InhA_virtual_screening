# HPC Next Steps

Upload this archive to the HPC:

```bash
md_top2_control_clean2aq8.zip
```

Then unpack it:

```bash
unzip md_top2_control_clean2aq8.zip
cd md_top2_control_clean2aq8
```

## Check available tools

Run:

```bash
module avail gromacs
module avail amber
module avail openbabel
module avail acpype
module avail charmm
module avail cgenff
```

If the cluster uses Lmod, also try:

```bash
module spider gromacs
module spider ambertools
module spider openbabel
module spider acpype
```

## Build one pilot system first

Start with:

```text
Lantadene_A
```

Then build:

```text
PT70_Control_Inhibitor
```

Only after those work, build:

```text
Rutin
```

## Required files before running the SLURM template

For each compound, you need a prepared simulation directory containing:

```text
complex.gro
topol.top
index.ndx
```

The included SLURM script expects those names.

Use the per-compound workspaces under:

```text
systems/
```

Detailed build guidance is in:

```text
HPC_BUILD_GUIDE.md
```

## Suggested force-field route

Use one consistent route for all three systems.

Recommended practical choices:

```text
Protein: CHARMM36m or AMBER ff14SB/ff19SB
Water: TIP3P
Ligands: CGenFF for CHARMM route, or GAFF/OpenFF for AMBER route
Ions: 0.15 M NaCl
Temperature: 300 K
Pressure: 1 bar
```

Do not mix unrelated force fields unless your workflow explicitly supports it.

## After the 10 ns pilot

Use:

```text
POST_RUN_ANALYSIS.md
```

to check RMSD, ligand retention, hydrogen bonding, and protein-ligand distance before committing to longer production MD.

## Key warning

The ligand PDBQT files are docking inputs/outputs. They are useful for poses, but they are not MD topology files. Use SDF/SMILES for ligand chemistry, then verify the generated ligand structure still matches the docked pose.

