ML4GutHealth
A Data-Driven Universal Gut Microbiome Health Assessment: A Machine Learning Framework trained on Large Metagenomic Data
Bablu Kumar1,2, Erika Lorusso2,3, Bruno Fosso2* and Graziano Pesole2,3*
 
1Department of Oncology and Hematology-Oncology, Università degli Studi di Milano, Milan, Italy.
2Department of Biosciences, Biotechnology and Environment, University of Bari A. Moro, Bari, Italy.
3National Research Council, Institute of Biomembranes, Bioenergetics and Molecular Biotechnologies, Bari, Italy.
 
*CORRESPONDENCE
Graziano Pesole: graziano.pesole@uniba.it; graziano.pesole@cnr.it
Bruno Fosso: bruno.fosso@uniba.it

Code and analysis scripts for a machine learning classifier that distinguishes healthy from non-healthy individuals using gut microbiome species-level profiles, trained on pooled public stool metagenomic data spanning multiple disease conditions.
Overview
This repository contains the processing pipeline, model training code, and supplementary analyses associated with the manuscript. The classifier was trained on 7,452 publicly available stool shotgun metagenomes from 32 independent studies, covering 12 disease conditions and healthy controls, and externally validated on 642 samples from 6 independent studies, including disease types not present in the training data.
Repository Structure
ML4GutHealth/
preprocessing/ computational_pipeline.sh 
diversity_analysis/alpha and betadivesity  alpha_diversity_calculation.R
‘- batch effects  analysis 
├── model_training/ all the model codes for the training 

Name	Last commit message	Last commit date
..
LR_ElasticNet_training.py
Add Logistic Regression training script with ElasticNet
1 minute ago
rf_model_training.py
Update comments for clarity on input filename
svm_model_trianing.py

├── sensitivity_analysis/
│   ├── normalization_strategy/ # Raw counts, TSS, CLR comparison
│   └── clr_pseudocount/        # CLR pseudocount sensitivity (1.0, 0.5, 0.1)
├── figures/                    # Scripts used to generate manuscript and supplementary figures
└── supplementary_tables/       # Supporting data tables

Data table: 
training_test_metadata_for_7452_32_studies.tsv
validation_metadata_642_6.csv
Manuscript figures 


Data
The species-level abundance matrices, associated sample metadata, and PFI-selected feature lists used in this repository are described in the manuscript's Data Availability section. Raw sequencing data were obtained from public repositories (SRA/ENA) under the accessions listed in Supplementary Table S2 (training) and Supplementary Table S7 (validation).
Pipeline Summary
Data retrieval: Raw reads downloaded using the SRA Toolkit.
Quality control: Fastp (v0.24.0), paired-end adapter detection, minimum quality score 20, minimum read length 45 bp.
Host-read removal: Bowtie2 (v2.5.4) against the GRCh38 no-alt reference genome.
Taxonomic classification: Kraken2 (v2.1.3) against the UHGG human-gut reference database (v2.0.2), refined with Bracken (v2.9) at species level.
Batch correction: MMUPHin, using BioProject ID as the batch variable.
Normalization: Centered log-ratio (CLR) transformation.
Model training: Four algorithms (SVM-Linear, SVM-RBF, LR-ElasticNet, Random Forest), each combined with three feature-selection strategies, trained and evaluated with stratified train/test splitting and hyperparameter tuning.
Feature selection: Permutation feature importance (PFI) applied to the best-performing model (SVM-RBF).
Interpretability: MaAsLin2 used to statistically confirm associations between PFI-selected species and health status.
Requirements
Python 3.10, conda environment (see environment.yml)
R (version specified in renv.lock), with packages: MMUPHin, MaAsLin2, iNEXT, ggplot2, patchwork
Key Python packages: scikit-learn, scikit-bio, pandas, numpy
Usage
Each subdirectory contains scripts intended to be run in the order listed in its own README or as numbered filenames. Input paths are set at the top of each script and should be adjusted to match the local data location.
Citation
If you use this code, please cite the associated manuscript (citation to be added upon publication).
Contact
For questions regarding this repository, please open an issue or contact the corresponding author.

