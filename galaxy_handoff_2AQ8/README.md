# Galaxy MD Handoff: Original 2AQ8 Compound Set

This package contains the inputs and results for the **original six-compound 2AQ8 docking set**. It is an upload/preparation bundle, not a completed molecular-dynamics system. Ligand parameters, protonation, solvent, ions, and simulation settings still need to be prepared and validated on the chosen Galaxy instance.

## Compounds and existing docking results

The six compounds are Lantadene A, Rutin, Guaijaverin, Quercitrin, 7-Methyljuglone, and the project's PT70 control. `results/final_docking_summary.csv` preserves the original AutoDock Vina scores and docking ligand efficiencies. `results/2aq8_pose_contact_summary.csv` and `results/2aq8_close_contact_pairs.csv` contain the original-pose geometric contact analysis.

## Structure files

- `structures/2AQ8_full_crystal.pdb`: deposited 2AQ8 chain-A crystal structure, including the bound NAD cofactor (PDB residue name `NAI`, chain A residue 270) and crystallographic waters/other records.
- `structures/2aq8_chainA_clean.pdb`: cleaned chain-A protein structure used as an input during project receptor preparation; protein-only, not a complete MD system.
- `structures/2aq8_prepared.pdb`: prepared receptor structure already present in the project.
- `structures/2aq8_receptor_docking_only.pdbqt`: receptor representation for docking only. Do not use PDBQT as the production MD topology.
- `structures/TCU_ideal.sdf`: RCSB ideal component structure for TCU from the supplemental 2X22 redocking check. This is NOT the project's original PT70 control structure and is included only as a reference; do not substitute it for the original control.

## Ligands

- `ligands/original_inputs/`: original PDBQT ligand inputs used for the prior 2AQ8 docking.
- `ligands/docked_poses/`: top-ranked Vina poses from those original 2AQ8 runs.
- `results/original_2AQ8_ligand_smiles.csv`: exact SMILES read from the original PDBQT input remarks, including the original project control structure.

PDBQT files retain docking coordinates and atom types, but are not parameterized MD ligand topologies. Use the corresponding original ligand structures/SMILES from the project's candidate data to create force-field-compatible ligand files, then map/check the prepared ligand against its docked pose. Explicitly verify stereochemistry, bond orders, protonation/tautomer state, and atom mapping before minimization.

## Critical MD preparation decisions

1. Decide whether NAD is retained. The original crystal includes NAD; the docking receptor used for the existing 2AQ8 screening removed residue `A:270` (NAI/NAD) and residue `A:447`. These deletions were made for docking and must not be inherited automatically by an MD system. For an enzyme complex, retaining and parameterizing NAD is likely the biologically relevant starting hypothesis, but confirm against the research question and force-field support.
2. Select and record protein and ligand force fields. Confirm the Galaxy workflow can parameterize each phytochemical and the PT70 control; do not proceed if any ligand lacks validated parameters.
3. Decide protonation states, termini, missing atoms/residues, crystallographic waters, ions, solvent model, and salt concentration consistently across systems.
4. Start with one original candidate and the original PT70 control as a preparation pilot. Minimize and equilibrate; check for severe clashes, energy/temperature/pressure stability, and ligand retention before expanding to all six.
5. Only after the pilot passes, run independent replicas with recorded random seeds and sufficient trajectory analysis/storage allocation.

## Suggested Galaxy workflow

1. Upload/unpack this archive in a Galaxy history.
2. Use `2AQ8_full_crystal.pdb` as the structural reference; inspect chain A and NAD. Compare to `2aq8_chainA_clean.pdb` only to understand prior docking preparation.
3. Prepare protein and decide NAD handling explicitly. Save the chosen prepared complex and parameterization inputs as new Galaxy datasets; leave these source files unchanged.
4. Prepare each ligand from verified chemical structures (not by treating PDBQT as a force-field topology) and place it into the chosen complex using the original docked pose coordinates after atom mapping.
5. Run a single-complex minimization/equilibration pilot, then production replicas only after checking system stability and available Galaxy quotas.
6. Analyze ligand/protein RMSD, pocket residue RMSF, contact persistence, and polar interaction occupancy. Treat these as simulation observations, not proof of inhibition.

## Reproducibility and scope

The prior Vina search used center `(2.2, -32.2, 14.0) A`, box size `(22, 22, 22) A`, exhaustiveness 12, and five modes. Existing original results are not overwritten by this Galaxy handoff. The separate `validation/exploratory_2X22/` panel is intentionally not included in this package.
