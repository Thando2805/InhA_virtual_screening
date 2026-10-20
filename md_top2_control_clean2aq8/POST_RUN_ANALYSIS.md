# Post-Run Analysis

Run these checks after each 10 ns pilot completes.

## Prepare a Centered Trajectory

```bash
gmx trjconv -s md_10ns_test.tpr -f md_10ns_test.xtc -o md_10ns_test_centered.xtc -pbc mol -center
```

Choose the protein or protein-ligand group for centering, then output the whole system or protein-ligand group as needed.

## Protein Backbone RMSD

```bash
gmx rms -s md_10ns_test.tpr -f md_10ns_test_centered.xtc -o rmsd_backbone.xvg
```

Select backbone for fitting and backbone for RMSD.

## Ligand RMSD

```bash
gmx rms -s md_10ns_test.tpr -f md_10ns_test_centered.xtc -n index.ndx -o rmsd_ligand.xvg
```

Use a ligand group if present in `index.ndx`.

## Radius of Gyration

```bash
gmx gyrate -s md_10ns_test.tpr -f md_10ns_test_centered.xtc -o gyrate.xvg
```

## Hydrogen Bonds

```bash
gmx hbond -s md_10ns_test.tpr -f md_10ns_test_centered.xtc -n index.ndx -num hbonds_protein_ligand.xvg
```

Use protein and ligand groups.

## Minimum Protein-Ligand Distance

```bash
gmx mindist -s md_10ns_test.tpr -f md_10ns_test_centered.xtc -n index.ndx -od mindist_protein_ligand.xvg
```

Use protein and ligand groups.

## Pass/Fail Signals for the Pilot

Continue to longer MD only if:

- Energy minimization converges without severe clashes.
- Temperature and pressure equilibrate normally.
- Protein backbone RMSD stabilizes instead of drifting continuously.
- Ligand remains in the intended InhA binding pocket.
- Protein-ligand distance and interaction checks do not show early dissociation.

