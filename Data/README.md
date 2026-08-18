Data Content
This folder contains the following files:

1. `raw_otu_count_table.csv.bz2`: feature table containing species raw relative abundance counts;
2. `training_test_metadata_for_7452_32_studies.tsv`: metadata for training and testing datasets;
3. `validation_metadata_642_6.csv`: metadata for the validation datasets;
4. `without_corrected_otu_7452_species.csv.bz2`: species-level OTU feature table for the 7452-sample training set, without batch correction applied;
5. `corrected_training_7452.csv.bz2`: species-level OTU feature table for the 7452-sample training set, after MMUPHin batch correction;
6. `without_corrected_validation_OTU_642_samples.csv.bz2`: species-level OTU feature table for the 642-sample validation set, without batch correction applied;
7. `corrected_corrected_validation_OTU_642_samples.csv.bz2`: species-level OTU feature table for the 642-sample validation set, after MMUPHin batch correction;
8. `batch_effect_analysis.ipynb`: notebook comparing corrected vs. uncorrected feature tables to assess batch effects across studies.
