# Research Progress Report

**Project:** Machine-learning-assisted virtual screening of Zimbabwean phytochemicals against *Mycobacterium tuberculosis* InhA  
**Date:** 7 October 2026  
**Prepared for:** Supervisor / Professor  
**Student:** [Add name]  
**Programme / course:** [Add programme or course]

## 1. Project Aim

This project is investigating whether selected phytochemicals associated with Zimbabwean medicinal plants may be worth prioritising for follow-up as potential inhibitors of enoyl-acyl carrier protein reductase (InhA), an enzyme from *M. tuberculosis*. The current work uses ChEMBL bioactivity records to train a molecular-structure classifier, applies that model to a curated phytochemical panel, and uses molecular docking to further rank selected candidates. The findings are computational predictions and have not yet been experimentally validated.

## 2. Work Completed

1. **Bioactivity data collection:** A raw activity table was assembled from ChEMBL records for InhA (target accession P9WGR1 / ChEMBL target CHEMBL1849). The saved table contains 457 IC50 records.
2. **Dataset processing:** The records were audited and restricted to positive, exact IC50 measurements reported in nM. This retained 402 activity records; 55 records were excluded (52 reported as greater-than values, 2 as less-than values, and 1 exact value reported in ug/mL). Duplicate canonical SMILES were consolidated by median pIC50, then labelled active at pIC50 >= 6 (1 micromolar or lower). The regenerated dataset contains 345 unique valid structures: 129 active and 216 inactive.
3. **Molecular descriptors and model training:** RDKit Morgan (ECFP4) fingerprints with 2,048 bits were generated. A class-weighted Random Forest was trained and evaluated on a stratified 80/20 split. XGBoost was unavailable in the training environment and was not included in this run.
4. **Phytochemical library assembly:** A structure-validated library of 44 entries was saved, representing compounds associated with 21 plant species plus a synthetic control entry.
5. **Candidate prediction and prioritisation:** The 25-compound panel was rescored with the newly trained Random Forest. PT70 was the only entry above the 0.50 classification threshold; none of the natural phytochemicals met that cutoff.
6. **Receptor and ligand preparation:** Project files include the InhA structure 2AQ8, a cleaned chain-A structure, prepared receptor files, and PDBQT ligand inputs. RDKit/Meeko preparation and AutoDock Vina docking scripts are present.
7. **Docking:** Six compounds selected using the earlier model, including the PT70 control, have saved docked poses and a ranked affinity summary.

## 3. Current Results

The ML values below are model scores from the current Random Forest, not calibrated probabilities of biological activity. Docking results are shown separately because the six docking runs were performed on compounds selected using the earlier model, not the current ranking. Vina scores are in kcal/mol; more negative scores rank more favourably within that docking run, but are not experimental binding free energies. Ligand efficiency is calculated as negative Vina score divided by heavy-atom count.

| Current ML rank | Compound | Plant source / role | Model score | Threshold result |
|---|---|---|---:|---|
| 1 | PT70 (control) | Synthetic control | 0.7875 | Above 0.50 |
| 2 | Rutin | *Ficus ingens* | 0.3585 | Below 0.50 |
| 3 | Guaijaverin | *Psidium guajava* | 0.3524 | Below 0.50 |
| 4 | Amentoflavone | *Searsia pyroides* | 0.3460 | Below 0.50 |
| 5 | Kaempferol | *Pellaea* sp. | 0.3311 | Below 0.50 |
| 6 | Quercitrin | *Piliostigma thonningii* | 0.3268 | Below 0.50 |

| Original 2AQ8 docking compound | Vina score (kcal/mol) | Heavy atoms | Ligand efficiency (kcal/mol/heavy atom) |
|---|---:|---:|---:|
| Lantadene A | -10.142 | 33 | 0.3073 |
| Rutin | -9.249 | 43 | 0.2151 |
| Guaijaverin | -8.908 | 30 | 0.2969 |
| Quercitrin | -8.762 | 32 | 0.2738 |
| 7-Methyljuglone | -7.999 | 14 | 0.5714 |
| PT70 (control entry) | -7.514 | 20 | 0.3757 |

Among the original 2AQ8 panel, Lantadene A had the most favourable docking score, while 7-Methyljuglone had the highest score per heavy atom. These docking measurements do not establish biological activity.

The top-pose contact screen found 9-18 protein residues with at least one ligand heavy-atom pair within 4.0 A per compound. Gly14 and Gly96 were within 4.0 A for all six compounds; Ser94 was close to five, and Phe41 to five. The number of atom pairs within 4.0 A ranged from 48 for 7-Methyljuglone to 98 for Rutin. N/O/S atom pairs within 3.5 A ranged from 2 for the PT70 control to 15 for Guaijaverin. These are geometric proximity counts only: they are not hydrogen-bond assignments and need visual/protonation/angle review. They also should not be compared as a direct measure of binding strength across different-sized ligands.

| Original 2AQ8 compound | Residues within 4.0 A | Heavy-atom pairs <=4.0 A | Close N/O/S pairs <=3.5 A |
|---|---:|---:|---:|
| Lantadene A | 17 | 67 | 5 |
| Rutin | 18 | 98 | 8 |
| Guaijaverin | 13 | 86 | 15 |
| Quercitrin | 17 | 58 | 9 |
| 7-Methyljuglone | 10 | 48 | 3 |
| PT70 control | 9 | 57 | 2 |

The newly trained Random Forest achieved ROC-AUC 0.8927, accuracy 0.8261, F1 0.7500, and MCC 0.6232 on a stratified random held-out test set of 69 compounds (276 compounds used for training; random state 42; decision threshold 0.50). Its confusion matrix was 39 true negatives, 4 false positives, 8 false negatives, and 18 true positives.

A separate Bemis-Murcko scaffold-disjoint holdout evaluation achieved ROC-AUC 0.9751, accuracy 0.9130, F1 0.8125, and MCC 0.7587 on 69 test compounds (276 training compounds). The split contained 116 training scaffolds and 30 test scaffolds, with zero scaffold overlap. It was selected from up to 1,000 deterministic group-shuffle candidates based on test-set size near 20% and presence of both activity classes, not model score.

Repeated scaffold-aware validation was then run using five repeats of five-fold StratifiedGroupKFold (25 folds total). Mean +/- standard deviation across folds was ROC-AUC 0.8619 +/- 0.1175, accuracy 0.8383 +/- 0.0941, F1 0.7525 +/- 0.1411, and MCC 0.6577 +/- 0.1969. Fold ranges were 0.5650-0.9791 for ROC-AUC, 0.5507-0.9565 for accuracy, 0.4151-0.9388 for F1, and 0.0525-0.9055 for MCC. All 25 folds had zero train/test scaffold overlap. The variability indicates that results depend materially on which scaffolds are held out; the earlier single scaffold holdout score should not be read as a stable expected performance estimate.

## 4. Interpretation and Limitations

- **No biological confirmation yet:** Neither the classifier nor docking establishes enzyme inhibition, antibacterial activity, selectivity, safety, or therapeutic potential. Laboratory testing is required.
- **Training results need further validation:** Metrics are now saved in `models/training_metrics.csv`, but they come from a single random split, not cross-validation or an external test set. Similar compounds may occur in both partitions.
- **Scaffold validation has substantial split variability:** Repeated five-fold scaffold-aware evaluation had zero train/test scaffold overlap in all folds, but ROC-AUC ranged from 0.5650 to 0.9791. The mean and standard deviation are more informative than the earlier favorable single split, and further data/model review is warranted before treating the classifier as robust.
- **Activity filtering is deliberately conservative:** Censored measurements and non-nM records are excluded rather than converted or interpreted. Review excluded records and assay context before deciding whether a justified unit conversion or censored-data method should recover any of them.
- **No natural candidate meets the current cutoff:** The top natural model score is 0.3585 for Rutin, below the project's 0.50 active threshold. These compounds should not be described as ML-predicted active hits; any further docking should be justified as exploratory ranking.
- **Current model panel differs from the original docked set:** The refreshed model's top natural candidates do not match the original six 2AQ8 docked compounds. The original docking set and its follow-up contact analysis are the focus of this report; scores from any separate receptor/panel run are not combined with them.
- **Control identity should be checked:** Different PT70 SMILES appear in the project scripts. Confirm the intended control structure and use the same verified structure for prediction, ligand preparation, and docking.
- **Docking protocol scope:** The original compounds were docked using the 2AQ8 receptor after the preparation choices recorded in `convert_receptor.py`. Confirm that the removed residues and cofactor state are appropriate before interpreting the poses biologically.
- **Candidate scope differs from the full library:** The saved plant library contains 44 entries, but the current ML screening script defines a 25-compound panel. The selection rationale and omitted compounds should be documented before describing the results as a complete screen of the library.
- **Input structures require verification:** Confirm compound identities, stereochemistry, and plant-source attribution against suitable chemical and phytochemical references before downstream experimental selection.

## 5. Next Steps

1. Visually inspect the six original 2AQ8 poses and verify the listed close polar geometries, including donor/acceptor identity and interaction angles.
2. Confirm the biological cofactor state intended for InhA docking and document the 2AQ8 receptor-preparation choices.
3. Investigate scaffold-fold model variability and preserve an external test set if additional activity data become available.
4. Plan in-vitro InhA testing of the original candidate set with positive, negative, and vehicle controls before making biological claims.

## 6. Project Files Supporting This Report

- `data/raw/inhA_chembl_ic50_raw.csv` — downloaded activity records.
- `data/processed/inhA_processed_ml.csv` — processed activity labels and fingerprints.
- `data/processed/inhA_dataset_audit.csv` — preprocessing retention and class-count audit.
- `data/raw/zimbabwe_phytochemicals.csv` — curated structure library.
- `data/processed/zimbabwean_phytochemical_hits.csv` — ML scores and candidate ranking.
- `models/best_inhA_model.joblib` — saved classifier.
- `models/training_metrics.csv` — held-out model metrics and confusion matrix.
- `models/scaffold_validation_metrics.csv` — scaffold-disjoint evaluation metrics.
- `models/scaffold_validation_split.csv` — molecule-level scaffold split assignments.
- `models/repeated_scaffold_cv_folds.csv` — per-fold repeated scaffold-aware cross-validation results.
- `models/repeated_scaffold_cv_summary.csv` — mean, standard deviation, and range across 25 folds.
- `final_docking_summary.csv` — ranked Vina results.
- `scripts/11_analyze_2aq8_pose_contacts.py` — heavy-atom proximity analysis for the original six 2AQ8 docking poses.
- `validation/original_2AQ8_contacts/2aq8_pose_contact_summary.csv` — per-compound residue-contact counts and closest residues.
- `validation/original_2AQ8_contacts/2aq8_close_contact_pairs.csv` — atom-pair distances for contacts within 4.0 A.
- Supplemental 2X22 validation and exploratory-panel files are retained under `validation/2X22_TCU_redocking/` and `validation/exploratory_2X22/`; they are separate from the primary 2AQ8 results summarized above.
- `scripts/01_fetch_chembl_data.py` through `scripts/05_fetch_dataset.py`, `run_vina.py`, and `convert_receptor.py` — workflow scripts.

**Overall status:** Activity preprocessing has been corrected to use exact IC50 values in nM, the Random Forest has been retrained and repeatedly evaluated with scaffold-aware splits, and the original six 2AQ8 docking poses have been analyzed for nearby protein contacts. The primary docking results remain the original 2AQ8 compounds; separate 2X22 exploratory work is retained outside this primary result set. No experimental confirmation has been performed.