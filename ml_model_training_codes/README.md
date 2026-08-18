
## Scripts

This folder contains the model training and evaluations script:

### `LR_ElasticNet_training.py`
Trains a Logistic Regression (ElasticNet) model on species-level, CLR-transformed abundance data, comparing 4 feature-selection strategies (All Features, RFECV, Lasso, PFI). Pass the prevalence-filtered input file as a command-line argument:
`python LR_ElasticNet_training.py <input_file.tsv>`

### `rf_model_training.py`
Trains a Random Forest classifier on the same CLR-transformed data, with the same 4 feature-selection strategies. Pass the input file as an argument:
`python rf_model_training.py clr_otu_metadata_merged_prevalence_bacterial_20pct.tsv`

### `svm_model_training.py`
Trains SVM models (both RBF and Linear kernels), comparing the same 4 feature-selection strategies. Pass the input file as an argument:
`python svm_model_training.py <input_file.tsv>`

### `feature_importance_svm_rbf_selective.py`
Runs after training — loads the pre-trained models and computes model-specific feature importance.

### `svm_rbf_pfi_normalization.py`
Checks whether normalization strategy (Raw, TSS, CLR with different pseudocounts, multiplicative zero-replacement) affects the performance of the best model (SVM-RBF, PFI features). Compares 6 configurations using fixed hyperparameters, with output saved to a single folder.

### `shap_on_models.py`
Runs a SHAP-based evaluation on a trained SVC model. Pass the model path as an argument:
`python shap_on_models.py -m <model_path>`

### `model_evaluation_full_pipeline.py`
Runs after training — evaluates all 4 trained models (RF, SVM-RBF, SVM-Linear, LogReg-ElasticNet) on the test set and the independent validation dataset (metrics, ROC/PR curves, confusion matrices, feature importance). Pass the prevalence value as an argument:
`python model_evaluation_full_pipeline.py 20pct`


### Note: All code was run as an HTCondor job on the ReCaS-Bari HPC cluster, with 32 CPU cores and 272 GB of RAM allocated per job.
