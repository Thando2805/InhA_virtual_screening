# PT70_Control_Inhibitor MD System Workspace

Use this directory to build the GROMACS-ready complex for the PT70 control run.

Expected files before submitting the SLURM job:

```text
complex.gro
topol.top
index.ndx
```

Copy or generate any force-field-specific ligand files here as needed. Keep the original docked pose in:

```text
../../ligands/PT70_Control_Inhibitor/PT70_Control_Inhibitor_docked.sdf
```

