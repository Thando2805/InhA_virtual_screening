# Rutin MD System Workspace

Use this directory to build the GROMACS-ready complex for the Rutin run after the pilot workflow is stable.

Expected files before submitting the SLURM job:

```text
complex.gro
topol.top
index.ndx
```

Copy or generate any force-field-specific ligand files here as needed. Keep the original docked pose in:

```text
../../ligands/Rutin/Rutin_docked.sdf
```

