# InhA Virtual Screening

Computational virtual screening pipeline for *Mycobacterium tuberculosis* Enoyl-ACP reductase (InhA), targeting novel phytochemical inhibitors from Zimbabwean medicinal plants.

## Overview

- **Target:** InhA (enoyl-acyl carrier protein reductase), PDB ID: 2AQ8
- **Approach:** ML-guided ligand prioritization + AutoDock Vina docking + ADMET/PKCSM prediction
- **Data Source:** Zimbabwean phytochemicals from ChEMBL / literature curation

## Repository Structure

```
inha_virtual_screening/
├── scripts/          # Python pipeline scripts (fetch, process, train, dock, analyze)
├── bin/              # Vina executables (excluded by .gitignore)
├── data/             # Large datasets (excluded by .gitignore)
├── structures/       # PDB/PDBQT structures (excluded by .gitignore)
├── models/           # Trained ML models (excluded by .gitignore)
├── ml_selected_pdbqt/# Ligand PDBQT files from ML selection
├── docking_results/  # Vina output poses (excluded by .gitignore)
├── project/          # Project-specific outputs (figures, tables, interactions)
├── galaxy_handoff_2AQ8/  # Galaxy server workflow files
├── md_hpc_start/     # GROMACS MD starting configs
├── md_top2_control_clean2aq8/  # MD top-2 control trajectories
└── validation/       # Validation results
```

## Pipeline

1. **Data Fetch** (`scripts/01_fetch_chembl_data.py`, `05_fetch_dataset.py`)
2. **Process & Label** (`scripts/02_process_and_label_data.py`, `04_virtual_screening.py`)
3. **Train ML Models** (`scripts/03_train_ml_models.py`, `06_scaffold_validation.py`, `07_repeated_scaffold_cv.py`)
4. **Validate Docking** (`scripts/08_validate_docking.py`, `09_exploratory_docking_2x22.py`, `10_exploratory_docking_2aq8.py`, `11_analyze_2aq8_pose_contacts.py`)
5. **Docking** (`run_vina.py`) — AutoDock Vina against 2AQ8
6. **Receptor Prep** (`convert_receptor.py`) — mk_prepare_receptor
7. **Analysis** (`finish_analysis.py`) — consolidated results + analysis

## Key Files

| File | Purpose |
| ---- | ------- |
| `zimbabwean_phytochemical_hits.csv` | Initial compound list (ML labels) |
| `final_docking_summary.csv` | Docking scores + ligand efficiency |
| `run_vina.py` | AutoDock Vina wrapper |
| `convert_receptor.py` | Receptor preparation |
| `finish_analysis.py` | Post-docking analysis |

## Environment

- Python 3.9+ with `venv` (`.venv/` excluded from Git)
- AutoDock Vina (`bin/vina.exe`)
- GROMACS (MD workflows)
- Required: `mk_prepare_receptor.exe`, `vina.exe` in `bin/`

## Prerequisites

1. Clone the repo
2. Create env: `python -m venv .venv && .venv\Scripts\activate`
3. Install: `pip install -r requirements.txt` (create if missing)
4. Place `vina.exe` and `mk_prepare_receptor.exe` in `bin/`
5. Place `2aq8_chainA_clean.pdb` in `data/structures/` (or adjust paths in `convert_receptor.py`)

## License

Proprietary research — do not redistribute.
