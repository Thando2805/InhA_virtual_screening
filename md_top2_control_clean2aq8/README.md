# MD Package: Top 2 Docking Hits + Control with Clean 2AQ8

This folder contains the selected systems for the molecular dynamics stage.

## Selection rule

The compounds were selected from `final_docking_summary.csv` by Vina binding affinity.

Selected systems:

| Role | Compound | Vina affinity, kcal/mol | Heavy atoms | Ligand efficiency |
| --- | --- | ---: | ---: | ---: |
| Top hit 1 | Lantadene_A | -10.142 | 33 | 0.3073 |
| Top hit 2 | Rutin | -9.249 | 43 | 0.2151 |
| Control | PT70_Control_Inhibitor | -7.514 | 20 | 0.3757 |

## Protein used

The protein file for this package is:

```text
protein/2aq8_chainA_clean.pdb
```

This is the clean 2AQ8 chain-A protein file requested for the MD stage.

## Ligand folders

Each selected ligand folder contains:

- Original PDBQT input used for docking.
- Top docked PDBQT pose.
- Converted docked SDF pose.

Ligand folders:

```text
ligands/Lantadene_A/
ligands/Rutin/
ligands/PT70_Control_Inhibitor/
```

## Important warning before MD

PDBQT files are docking files, not molecular dynamics topologies. For MD, generate force-field-compatible ligand parameters from verified ligand chemistry/SMILES/SDF, then align/check the prepared ligand against the docked pose.

Do not run production MD until each system has:

```text
complex.gro
topol.top
index.ndx
```

## Recommended order on the HPC

1. Build `Lantadene_A` complex first.
2. Build `PT70_Control_Inhibitor` control second.
3. Run minimization, NVT, NPT, and a short 10 ns pilot.
4. Check stability and ligand retention.
5. Build and run `Rutin` after the pilot workflow is stable.

## Included templates

GROMACS MDP files:

```text
mdp/minim.mdp
mdp/nvt.mdp
mdp/npt.mdp
mdp/md_10ns_test.mdp
```

SLURM template:

```text
slurm/gromacs_pipeline_template.sbatch
```

Edit the partition name, module name, and resource requests before submitting on the HPC.

## Added build scaffolding

Per-system workspaces are provided in:

```text
systems/Lantadene_A/
systems/PT70_Control_Inhibitor/
systems/Rutin/
```

Use `HPC_BUILD_GUIDE.md` to generate `complex.gro`, `topol.top`, and `index.ndx` for each system before submitting the SLURM job.

After a 10 ns pilot completes, use `POST_RUN_ANALYSIS.md` for the first stability checks.

