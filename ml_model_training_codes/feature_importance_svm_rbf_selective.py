#! /lustrehome/babluuniba2022/miniconda3/envs/ml/bin/python

# ==============================================================================
# Feature Importance Extraction Pipeline
#
# Position in the overall workflow: run AFTER model training. Loads the
# pre-trained model files (.joblib) for all four algorithms (RF, SVM-Linear,
# SVM-RBF, LogReg-ElasticNet) across all four feature-selection strategies
# (All Features, RFECV, PFI, LASSO), and computes model-appropriate feature
# importance for each:
#   - RF: built-in impurity-based feature_importances_
#   - SVM-Linear / LogReg: model coefficients (coef_)
#   - SVM-RBF: permutation importance (F1-scored, 5 repeats), since RBF-SVM
#     has no coef_ attribute
#
# NOTE: In this run, only two calls are active:
#   - Test dataset:       SVM_RBF, PFI-selected features only
#   - Validation dataset: SVM_RBF, All Features only
# All other model / feature-selection combinations are commented out below
# but retained in the code for reference and future re-activation.
#
# Importance is computed on both the held-out test set and the independent
# validation dataset (n = 642, 6 studies), overall and per BioProject.
# ==============================================================================

import os
import sklearn
from sklearn.model_selection import train_test_split
import csv
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
from sklearn.metrics import (accuracy_score, precision_score, recall_score, f1_score, balanced_accuracy_score,
                             roc_auc_score, average_precision_score, precision_recall_curve, roc_curve, auc)
from sklearn.metrics import confusion_matrix
from sklearn.inspection import permutation_importance
import joblib
from sklearn.model_selection import learning_curve
from sklearn.metrics import classification_report

print("NumPy version:", np.__version__)
print("scikit-learn version:", sklearn.__version__)
print("Pandas version:", pd.__version__)

data_dir = "/lustrehome/babluuniba2022/bablu/condor/Data_Analysis/JUNE_2025/prevalence_based_threshold"
svm_model_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/svm_models_full"
rf_model_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/random_forest_models"
logreg_model_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/logreg_models"
lasso_feature_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/lasso_models"
rfecv_feature_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/svm_models_full"
pfi_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev"
output_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/2026_featrues_importance"
featrue_importance_dir = os.path.join(output_dir, "featrue_importance")
roc_data_dir = os.path.join(output_dir, "roc_data")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(roc_data_dir, exist_ok=True)
os.makedirs(featrue_importance_dir, exist_ok=True)

if len(sys.argv) < 2:
    print("Usage: python code.py <prevalence_value(e.g., '30pct')>")
    raise ValueError("Please provide a prevalence value (e.g., 'python validation.py 30pct')")
prevalence = sys.argv[1]

print(f"Loading models for prevalence: {prevalence}")
models_all_features = {
    # "RF": joblib.load(os.path.join(rf_model_dir, f'RF_BACTERIAL_bacterial_{prevalence}_ALL_FEATURES.joblib')),
    #"SVM_LIN": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_ALL_FEATURES_linear.joblib')),
    "SVM_RBF": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_ALL_FEATURES_rbf.joblib')),
    #"LogReg": joblib.load(os.path.join(logreg_model_dir, f'LOGREG_BACTERIAL_{prevalence}_ALL_FEATURES.joblib'))
}
models_pfi = {
    # "RF_PFI": joblib.load(os.path.join(rf_model_dir, f'RF_BACTERIAL_bacterial_{prevalence}_PERMUTATION_FEATURES.joblib')),
    #"SVM_LIN_PFI": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_PERMUTATION_FEATURES_linear.joblib')),
    "SVM_RBF_PFI": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_PERMUTATION_FEATURES_rbf.joblib')),
    #"LogReg_PFI": joblib.load(os.path.join(logreg_model_dir, f'LOGREG_BACTERIAL_{prevalence}_PERMUTATION_FEATURES.joblib'))
}
models_rfe = {
    # "RF_RFECV": joblib.load(os.path.join(rf_model_dir, f'RF_BACTERIAL_{prevalence}_RFECV_FEATURES.joblib')),
    #"SVM_LIN_RFECV": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_RFECV_FEATURES_linear.joblib')),
    "SVM_RBF_RFECV": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_RFECV_FEATURES_rbf.joblib')),
    #"LogReg_RFECV": joblib.load(os.path.join(logreg_model_dir, f'LOGREG_BACTERIAL_{prevalence}_RFECV_FEATURES.joblib'))
}
models_lasso = {
    # "RF_LASSO": joblib.load(os.path.join(rf_model_dir, f'RF_BACTERIAL_{prevalence}_LASSO_FEATURES.joblib')),
    #"SVM_LIN_LASSO": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_LASSO_FEATURES_linear.joblib')),
    "SVM_RBF_LASSO": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_LASSO_FEATURES_rbf.joblib')),
    #"LogReg_LASSO": joblib.load(os.path.join(logreg_model_dir, f'LOGREG_BACTERIAL_{prevalence}_LASSO_FEATURES.joblib'))
}

print(f"\nLoading test dataset for prevalence: {prevalence}")
test_dataset_path = os.path.join(data_dir, f'clr_otu_metadata_merged_prevalence_bacterial_{prevalence}.tsv')
df = pd.read_csv(test_dataset_path, sep='\t', low_memory=False)
df_copy = df.copy()
columns_to_drop = ["Run_ID", "BioProject", "BioSample", "Health_status", "Phenotype", "Full_Name", "Sex", "Age", "Location", "Sample", "BMI", "Fit", "Platform", "Author"]
df = df.drop(columns=[col for col in columns_to_drop if col in df.columns])
y = df['Class_Label']
X = df.drop(columns=['Class_Label'])
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
print(f"Dimension:\nTrain Rows:{X_train.shape[0]} | Features:{X_train.shape[1]}\nTest Rows:{X_test.shape[0]} | Features:{X_test.shape[1]}")
print("Overall phenotypes in the training and test dataset:\n", df_copy["Phenotype"].value_counts())

train_indices = X_train.index
test_indices = X_test.index
phenotypes_train = df_copy.loc[train_indices, 'Full_Name']
phenotypes_test = df_copy.loc[test_indices, 'Full_Name']
print(f"Number of unique phenotypes in training set: {phenotypes_train.nunique()}")
print(f"Number of unique phenotypes in test set: {phenotypes_test.nunique()}")
print("Phenotype counts in training set:\n", phenotypes_train.value_counts())
print("\nPhenotype counts in test set:\n", phenotypes_test.value_counts())

validation_dataset_path = "/lustrehome/babluuniba2022/bablu/condor/Data_Analysis/JUNE_2025/FEB_2026/clr_transformed_validation_dataset_test_642_6studies.csv"
dataset = pd.read_csv(validation_dataset_path, sep=',', low_memory=False)
print("\nDisease included in validation dataset:\n", dataset["Phenotype"].value_counts())
copy_df = dataset.copy()
validation_df = dataset.drop(columns=[col for col in columns_to_drop if col in dataset.columns])
y_validation = validation_df['Class_Label'].astype(int)
X_validation_all = validation_df.drop(columns=['Class_Label'])


rfe_df = pd.read_csv(os.path.join(rfecv_feature_dir, f'RFECV_SELECTED_FEATURES_WITH_COEFFICIENTS_{prevalence}.tsv'), sep='\t')
pfi = pd.read_csv(os.path.join(pfi_dir, f'PERMUTATION_IMPORTANCE_BACTERIAL_bacterial_{prevalence}.tsv'), sep='\t')
lasso_species = pd.read_csv(os.path.join(lasso_feature_dir, f'lasso_based_featrues_on_{prevalence}.tsv'), sep='\t')

selected_features = X_train.columns.tolist()
X_validation = X_validation_all[[col for col in selected_features if col in X_validation_all.columns]]

rfe_features = rfe_df["Feature_Name"].tolist()
X_validation_rfe = X_validation_all[[col for col in rfe_features if col in X_validation_all.columns]]
X_validation_rfe = X_validation_rfe.dropna()
y_validation_rfe = y_validation.loc[X_validation_rfe.index]

pfi_features = pfi[pfi['Importance_Mean'] > 0]["Feature"].tolist()
X_validation_pfi = X_validation_all[[col for col in pfi_features if col in X_validation_all.columns]]
y_validation_pfi = y_validation.loc[X_validation_pfi.index]

lasso_features = lasso_species["Feature"].tolist()
X_validation_lasso = X_validation_all[[col for col in lasso_features if col in X_validation_all.columns]]
y_validation_lasso = y_validation.loc[X_validation_lasso.index]


def evaluate_feature_importance(models, X, y, feature_selection_method, dataset_type, output_dir, bioproject_id=None):
    feature_importance_list = []
    for model_name, model in models.items():
        print(f"Processing model: {model_name}, type: {type(model)}")
        importance_df = None

        if len(np.unique(y)) < 2:
            print(f"Skipping feature importance for {model_name} in {feature_selection_method} (Dataset: {dataset_type}): only one class present")
            continue
        if 'LogReg' in model_name:
            print(f"Using coef_ for {model_name}")
            try:
                importance_df = pd.DataFrame({'Feature': X.columns, 'Importance': model.coef_.flatten(), 'Model': model_name,
                                              'Feature_Selection': feature_selection_method, 'Prevalence': prevalence, 'Dataset': dataset_type})
            except AttributeError as e:
                print(f"Error computing coef_ for {model_name}: {e}")
                continue
        elif 'SVM_RBF' in model_name:
            print(f"Using permutation_importance for {model_name}")
            try:
                perm_importance = permutation_importance(model, X, y, scoring='f1', n_repeats=5, random_state=42, n_jobs=-1)
                importance = perm_importance.importances_mean
                importance_df = pd.DataFrame({'Feature': X.columns, 'Importance': importance, 'Model': model_name,
                                              'Feature_Selection': feature_selection_method, 'Prevalence': prevalence, 'Dataset': dataset_type})
            except Exception as e:
                print(f"Error computing permutation_importance for {model_name}: {e}")
                continue
        elif 'SVM_LIN' in model_name:
            print(f"Using coef_ for {model_name}")
            try:
                importance = model.coef_.flatten()
                importance_df = pd.DataFrame({'Feature': X.columns, 'Importance': importance, 'Model': model_name, 'Feature_Selection': feature_selection_method,
                                              'Prevalence': prevalence, 'Dataset': dataset_type})
            except AttributeError as e:
                print(f"Error computing coef_ for {model_name}: {e}")
                continue
        elif 'RF' in model_name:
            print(f"Using feature_importances_ for {model_name}")
            try:
                importance = model.feature_importances_
                importance_df = pd.DataFrame({'Feature': X.columns, 'Importance': importance, 'Model': model_name, 'Feature_Selection': feature_selection_method,
                                              'Prevalence': prevalence, 'Dataset': dataset_type})
            except AttributeError as e:
                print(f"Error computing feature_importances_ for {model_name}: {e}")
                continue
        else:
            print(f"Skipping unknown model type: {model_name}")
            continue
        if bioproject_id:
            importance_df['BioProject'] = bioproject_id
        if importance_df is not None:
            output_file = os.path.join(output_dir, f'feature_importance_{model_name}_{feature_selection_method}_'
                f'{"bioproject_" + bioproject_id + "_" if bioproject_id else ""}{dataset_type}_{prevalence}.tsv')
            importance_df.to_csv(output_file, sep='\t', index=False)
            print(f"Feature importance saved to: {output_file}")
            feature_importance_list.append(importance_df)
        else:
            print(f"No feature importance generated for {model_name} in {feature_selection_method} (Dataset: {dataset_type})")
    return feature_importance_list


print("\nEvaluating feature importance on test dataset...")
X_test_rfe = X_test[[col for col in rfe_features if col in X_test.columns]]
X_test_pfi = X_test[[col for col in pfi_features if col in X_test.columns]]
X_test_lasso = X_test[[col for col in lasso_features if col in X_test.columns]]

all_feature_importance = []
# --- TEST DATASET: only SVM_RBF_PFI is active; all other calls kept for reference ---
# all_feature_importance.extend(evaluate_feature_importance(models_all_features, X_test, y_test, "All_Features", "test", featrue_importance_dir))
# all_feature_importance.extend(evaluate_feature_importance(models_rfe, X_test_rfe, y_test, "RFE", "test", featrue_importance_dir))
all_feature_importance.extend(evaluate_feature_importance(models_pfi, X_test_pfi, y_test, "PFI", "test", featrue_importance_dir))
# all_feature_importance.extend(evaluate_feature_importance(models_lasso, X_test_lasso, y_test, "LASSO", "test", featrue_importance_dir))

if all_feature_importance:
    pd.concat(all_feature_importance, ignore_index=True).to_csv(os.path.join(featrue_importance_dir, f'feature_importance_test_dataset_{prevalence}.tsv'), sep='\t', index=False)
    print(f"Combined feature importance saved to: {os.path.join(featrue_importance_dir, f'feature_importance_test_dataset_{prevalence}.tsv')}")

print("\nEvaluating feature importance on validation dataset...")
all_feature_importance = []
# --- VALIDATION DATASET: only SVM_RBF, All_Features is active; all other calls kept for reference ---
all_feature_importance.extend(evaluate_feature_importance(models_all_features, X_validation, y_validation, "All_Features", "validation", featrue_importance_dir))
# all_feature_importance.extend(evaluate_feature_importance(models_rfe, X_validation_rfe, y_validation_rfe, "RFE", "validation", featrue_importance_dir))
# all_feature_importance.extend(evaluate_feature_importance(models_pfi, X_validation_pfi, y_validation_pfi, "PFI", "validation", featrue_importance_dir))
# all_feature_importance.extend(evaluate_feature_importance(models_lasso, X_validation_lasso, y_validation_lasso, "LASSO", "validation", featrue_importance_dir))

if all_feature_importance:
    pd.concat(all_feature_importance, ignore_index=True).to_csv(os.path.join(featrue_importance_dir, f'feature_importance_validation_dataset_{prevalence}.tsv'), sep='\t', index=False)
    print(f"Combined feature importance saved to: {os.path.join(featrue_importance_dir, f'feature_importance_validation_dataset_{prevalence}.tsv')}")

# BioProject-based feature importance from test dataset
print("\nEvaluating feature importance on each BioProject in test datasets...")
test_grouped_data = df_copy.loc[test_indices].groupby('BioProject')
all_bioproject_feature_importance_test = []
for project_id, group in test_grouped_data:
    phenotype = group["Phenotype"].unique()
    phenotype_str = ",".join(phenotype)
    phenotype_counts = group["Phenotype"].value_counts().to_dict()
    phenotype_counts_str = f"{project_id}: {','.join([f'{key}({value})' for key, value in phenotype_counts.items()])}"

    test_df_bioproject = group.drop(columns=[col for col in columns_to_drop if col in group.columns])
    y_test_bioproject = test_df_bioproject['Class_Label'].astype(int)
    X_test_all = test_df_bioproject.drop(columns=['Class_Label'])

    X_test_bioproject = X_test_all[[col for col in X_train.columns if col in X_test_all.columns]]
    X_test_rfe = X_test_all[[col for col in rfe_features if col in X_test_all.columns]]
    X_test_rfe = X_test_rfe.dropna()
    X_test_pfi = X_test_all[[col for col in pfi_features if col in X_test_all.columns]]
    X_test_lasso = X_test_all[[col for col in lasso_features if col in X_test_all.columns]]

    y_test_rfe = y_test_bioproject.loc[X_test_rfe.index]
    y_test_pfi = y_test_bioproject.loc[X_test_pfi.index]
    y_test_lasso = y_test_bioproject.loc[X_test_lasso.index]

    print(f"BioProject {project_id} - y_test_bioproject unique values: {np.unique(y_test_bioproject)}")
    print(f"BioProject {project_id} - X_test_bioproject shape: {X_test_bioproject.shape}")
    print(f"BioProject {project_id} - X_test_rfe shape: {X_test_rfe.shape}")
    print(f"BioProject {project_id} - y_test_rfe shape: {y_test_rfe.shape}")

    # --- TEST DATASET (per BioProject): only SVM_RBF_PFI is active; all other calls kept for reference ---
    # all_bioproject_feature_importance_test.extend(evaluate_feature_importance(models_all_features, X_test_bioproject, y_test_bioproject, "All_Features", "bioproject", featrue_importance_dir, project_id))
    # all_bioproject_feature_importance_test.extend(evaluate_feature_importance(models_rfe, X_test_rfe, y_test_rfe, "RFE", "bioproject", featrue_importance_dir, project_id))
    all_bioproject_feature_importance_test.extend(evaluate_feature_importance(models_pfi, X_test_pfi, y_test_pfi, "PFI", "bioproject", featrue_importance_dir, project_id))
    # all_bioproject_feature_importance_test.extend(evaluate_feature_importance(models_lasso, X_test_lasso, y_test_lasso, "LASSO", "bioproject", featrue_importance_dir, project_id))

if all_bioproject_feature_importance_test:
    pd.concat(all_bioproject_feature_importance_test, ignore_index=True).to_csv(os.path.join(output_dir, f'feature_importance_bioproject_test_{prevalence}.tsv'), sep='\t', index=False)
    print(f"Combined test BioProject feature importance saved to: {os.path.join(output_dir, f'feature_importance_bioproject_test_{prevalence}.tsv')}")

# BioProject-based feature importance from validation dataset
print("\nEvaluating feature importance on each BioProject in validation datasets...")
grouped_data = dataset.groupby('BioProject')
all_bioproject_feature_importance = []
for project_id, group in grouped_data:
    phenotype = group["Phenotype"].unique()
    phenotype_str = ",".join(phenotype)
    phenotype_counts = group["Phenotype"].value_counts().to_dict()
    phenotype_counts_str = f"{project_id}\n{','.join([f'{key}({value})' for key, value in phenotype_counts.items()])}"
    print(f"\nProcessing BioProject: {project_id}\tPhenotype: {phenotype_str}")

    validation_df = group.drop(columns=[col for col in columns_to_drop if col in group.columns])
    y_validation = validation_df['Class_Label'].astype(int)
    X_validation_all = validation_df.drop(columns=['Class_Label'])

    X_validation = X_validation_all[[col for col in X_train.columns if col in X_validation_all.columns]]
    X_validation_rfe = X_validation_all[[col for col in rfe_features if col in X_validation_all.columns]]
    X_validation_rfe = X_validation_rfe.dropna()
    X_validation_pfi = X_validation_all[[col for col in pfi_features if col in X_validation_all.columns]]
    X_validation_lasso = X_validation_all[[col for col in lasso_features if col in X_validation_all.columns]]

    y_validation_rfe = y_validation.loc[X_validation_rfe.index]
    y_validation_pfi = y_validation.loc[X_validation_pfi.index]
    y_validation_lasso = y_validation.loc[X_validation_lasso.index]

    print(f"BioProject {project_id} - y_validation unique values: {np.unique(y_validation)}")
    print(f"BioProject {project_id} - X_validation shape: {X_validation.shape}")
    print(f"BioProject {project_id} - X_validation_rfe shape: {X_validation_rfe.shape}")
    print(f"BioProject {project_id} - y_validation_rfe shape: {y_validation_rfe.shape}")

    # --- VALIDATION DATASET (per BioProject): only SVM_RBF, All_Features is active; all other calls kept for reference ---
    all_bioproject_feature_importance.extend(evaluate_feature_importance(models_all_features, X_validation, y_validation, "All_Features", "bioproject", featrue_importance_dir, project_id))
    # all_bioproject_feature_importance.extend(evaluate_feature_importance(models_rfe, X_validation_rfe, y_validation_rfe, "RFE", "bioproject", featrue_importance_dir, project_id))
    # all_bioproject_feature_importance.extend(evaluate_feature_importance(models_pfi, X_validation_pfi, y_validation_pfi, "PFI", "bioproject", featrue_importance_dir, project_id))
    # all_bioproject_feature_importance.extend(evaluate_feature_importance(models_lasso, X_validation_lasso, y_validation_lasso, "LASSO", "bioproject", featrue_importance_dir, project_id))

if all_bioproject_feature_importance:
    pd.concat(all_bioproject_feature_importance, ignore_index=True).to_csv(os.path.join(output_dir, f'feature_importance_bioproject_{prevalence}.tsv'), sep='\t', index=False)
    print(f"Combined BioProject feature importance saved to: {os.path.join(output_dir, f'feature_importance_bioproject_{prevalence}.tsv')}")