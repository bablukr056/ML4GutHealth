#!/lustrehome/babluuniba2022/miniconda3/envs/ml/bin/python
"""
==============================================================================
Model Evaluation Pipeline - Prevalence-Filtered Feature Set (Run: 20% Prevalence)
==============================================================================

Position in the overall workflow
---------------------------------
This script is run AFTER model training and feature selection 

Purpose
-------
Performs the full downstream evaluation of the four base classifiers
(SVM-RBF, SVM-Linear, Random Forest, Logistic Regression / ElasticNet)
trained on bacterial species abundance clr transformed data 
threshold, across all four feature-selection strategies:
    - All Features
    - LASSO-selected features
    - Permutation Feature Importance (PFI)-selected features
    - Recursive Feature Elimination with Cross-Validation (RFECV)-selected features

For each feature-selection strategy and each model, the script computes:
    - Performance metrics (Accuracy, Balanced Accuracy, Precision, Recall,
      F1 Score, AUC, AUPR) on the held-out test split and on the
      independent validation cohort (n = 642, 6 studies)
    - ROC and Precision-Recall curves (overall, and per BioProject)
    - Confusion matrices (counts and percentages)
    - Per-sample classification outcomes and misclassification analysis
    - Per-class classification reports (precision / recall / F1 per class)
    - Permutation- / coefficient-based feature importance per model

This particular run evaluates the models trained on the 20% prevalence-
filtered feature table (prevalence argument = "20pct"), used to generate
the main-text and supplementary performance figures (Figures 2, 3, 4, 5;
Supplementary Figure 6).

Execution environment
----------------------
This script is intended to run as a non-interactive batch job on the HPC
cluster (Condor job submission) and NOT interactively, since it
loads multiple pre-trained models and processes the full dataset per
BioProject, which is memory- and time-intensive.

Usage:
    python evaluate_models.py <prevalence_value>
    e.g.
    python evaluate_models.py 20pct

Inputs
------
    - Pre-trained model files (.joblib) for RF, SVM-Linear, SVM-RBF, and
      LogReg-ElasticNet, for each feature-selection strategy (produced by
      the model-training scripts listed above)
    - Prevalence-filtered test dataset (CLR-transformed OTU table + metadata)
    - Independent validation dataset (n = 642, 6 studies)
    - LASSO / RFECV / PFI selected-feature lists for the given prevalence

Outputs
-------
All outputs are written under <output_dir>, including:
    - Performance_on_test_dataset_<prevalence>.tsv
    - Performance_on_validation_dataset_<prevalence>.tsv
    - classification_reports/          (per-model, per-feature-set, per-BioProject)
    - confusion_matrix_results/        (counts, plots, misclassifications)
    - roc_data/                        (ROC / PR curve data and facet plots)
    - featrue_importance/              (per-model feature importance tables)
    - model_evaluation_on_each_bioproject_test_datasets_<prevalence>.csv
    - model_evaluation_on_each_bioproject_validation_datasets_<prevalence>.csv

Note: the script is generic with respect to prevalence threshold and can be
re-run for any other value by changing the command-line argument; this
listing documents the 20% prevalence run specifically.
==============================================================================
"""

import os
import sys

import numpy as np
import pandas as pd
import joblib
import sklearn

from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    balanced_accuracy_score, roc_auc_score, average_precision_score,
    precision_recall_curve, roc_curve, auc, confusion_matrix,
    classification_report
)
from sklearn.inspection import permutation_importance

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

print("NumPy version:", np.__version__)
print("scikit-learn version:", sklearn.__version__)
print("Pandas version:", pd.__version__)


# ==============================================================================
# CONFIGURATION: paths, directories, prevalence argument
# ==============================================================================
data_dir = "/lustrehome/babluuniba2022/bablu/condor/Data_Analysis/JUNE_2025/prevalence_based_threshold"
svm_model_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/svm_models_full"
rf_model_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/random_forest_models"
logreg_model_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/logreg_models"
lasso_feature_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/lasso_models"
rfecv_feature_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/svm_models_full"
pfi_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev"
output_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/2026_perfromacne_20pct"

featrue_importance_dir = os.path.join(output_dir, "featrue_importance")
roc_data_dir = os.path.join(output_dir, "roc_data")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(roc_data_dir, exist_ok=True)
os.makedirs(featrue_importance_dir, exist_ok=True)

if len(sys.argv) < 2:
    print("Usage: python evaluate_models.py <prevalence_value> (e.g. '20pct')")
    raise ValueError("Please provide a prevalence value (e.g., 'python evaluate_models.py 20pct')")
prevalence = sys.argv[1]

columns_to_drop = [
    "Run_ID", "BioProject", "BioSample", "Health_status", "Phenotype", "Full_Name",
    "Sex", "Age", "Location", "Sample", "BMI", "Fit", "Platform", "Author"
]


# ==============================================================================
# LOAD PRE-TRAINED MODELS (one set per feature-selection strategy)
# ==============================================================================
print(f"Loading models for prevalence: {prevalence}")

models_all_features = {
    "RF": joblib.load(os.path.join(rf_model_dir, f'RF_BACTERIAL_bacterial_{prevalence}_ALL_FEATURES.joblib')),
    "SVM_LIN": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_ALL_FEATURES_linear.joblib')),
    "SVM_RBF": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_ALL_FEATURES_rbf.joblib')),
    "LogReg": joblib.load(os.path.join(logreg_model_dir, f'LOGREG_BACTERIAL_{prevalence}_ALL_FEATURES.joblib'))
}
models_pfi = {
    "RF_PFI": joblib.load(os.path.join(rf_model_dir, f'RF_BACTERIAL_bacterial_{prevalence}_PERMUTATION_FEATURES.joblib')),
    "SVM_LIN_PFI": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_PERMUTATION_FEATURES_linear.joblib')),
    "SVM_RBF_PFI": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_PERMUTATION_FEATURES_rbf.joblib')),
    "LogReg_PFI": joblib.load(os.path.join(logreg_model_dir, f'LOGREG_BACTERIAL_{prevalence}_PERMUTATION_FEATURES.joblib'))
}
models_rfe = {
    "RF_RFECV": joblib.load(os.path.join(rf_model_dir, f'RF_BACTERIAL_{prevalence}_RFECV_FEATURES.joblib')),
    "SVM_LIN_RFECV": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_RFECV_FEATURES_linear.joblib')),
    "SVM_RBF_RFECV": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_RFECV_FEATURES_rbf.joblib')),
    "LogReg_RFECV": joblib.load(os.path.join(logreg_model_dir, f'LOGREG_BACTERIAL_{prevalence}_RFECV_FEATURES.joblib'))
}
models_lasso = {
    "RF_LASSO": joblib.load(os.path.join(rf_model_dir, f'RF_BACTERIAL_{prevalence}_LASSO_FEATURES.joblib')),
    "SVM_LIN_LASSO": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_LASSO_FEATURES_linear.joblib')),
    "SVM_RBF_LASSO": joblib.load(os.path.join(svm_model_dir, f'SVM_BACTERIAL_{prevalence}_LASSO_FEATURES_rbf.joblib')),
    "LogReg_LASSO": joblib.load(os.path.join(logreg_model_dir, f'LOGREG_BACTERIAL_{prevalence}_LASSO_FEATURES.joblib'))
}


# ==============================================================================
# LOAD TEST DATASET AND REPRODUCE THE TRAIN/TEST SPLIT (seed=42, same as training)
# ==============================================================================
print(f"\nLoading test dataset for prevalence: {prevalence}")
test_dataset_path = os.path.join(data_dir, f'clr_otu_metadata_merged_prevalence_bacterial_{prevalence}.tsv')
df = pd.read_csv(test_dataset_path, sep='\t', low_memory=False)
df_copy = df.copy()

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


# ==============================================================================
# LOAD INDEPENDENT VALIDATION DATASET (n = 642, 6 studies)
# ==============================================================================
validation_dataset_path = "/lustrehome/babluuniba2022/bablu/condor/Data_Analysis/JUNE_2025/FEB_2026/clr_transformed_validation_dataset_test_642_6studies.csv"
dataset = pd.read_csv(validation_dataset_path, sep=',', low_memory=False)
print("\nDisease included in validation dataset:\n", dataset["Phenotype"].value_counts())

copy_df = dataset.copy()
validation_df = dataset.drop(columns=[col for col in columns_to_drop if col in dataset.columns])
y_validation = validation_df['Class_Label'].astype(int)
X_validation_all = validation_df.drop(columns=['Class_Label'])


# ==============================================================================
# LOAD FEATURE-SELECTION LISTS (RFECV / PFI / LASSO) AND SUBSET THE VALIDATION SET
# ==============================================================================
rfe_df = pd.read_csv(os.path.join(rfecv_feature_dir, f'RFECV_SELECTED_FEATURES_WITH_COEFFICIENTS_{prevalence}.tsv'), sep='\t')
pfi = pd.read_csv(os.path.join(pfi_dir, f'PERMUTATION_IMPORTANCE_BACTERIAL_bacterial_{prevalence}.tsv'), sep='\t')
lasso_species = pd.read_csv(os.path.join(lasso_feature_dir, f'lasso_based_featrues_on_{prevalence}.tsv'), sep='\t')

selected_features = X_train.columns.tolist()
X_validation = X_validation_all[[col for col in selected_features if col in X_validation_all.columns]]

rfe_features = rfe_df["Feature_Name"].tolist()
X_validation_rfe = X_validation_all[[col for col in rfe_features if col in X_validation_all.columns]].dropna()
y_validation_rfe = y_validation.loc[X_validation_rfe.index]

pfi_features = pfi[pfi['Importance_Mean'] > 0]["Feature"].tolist()
X_validation_pfi = X_validation_all[[col for col in pfi_features if col in X_validation_all.columns]]
y_validation_pfi = y_validation.loc[X_validation_pfi.index]

lasso_features = lasso_species["Feature"].tolist()
X_validation_lasso = X_validation_all[[col for col in lasso_features if col in X_validation_all.columns]]
y_validation_lasso = y_validation.loc[X_validation_lasso.index]


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def evaluate_model(model, X_test, y_test, y_prob=None):
    """Compute standard classification metrics for a single fitted model."""
    y_pred = model.predict(X_test)
    metrics = {
        'Accuracy': round(accuracy_score(y_test, y_pred), 4),
        'Balanced Accuracy': round(balanced_accuracy_score(y_test, y_pred), 4),
        'Precision': round(precision_score(y_test, y_pred, zero_division=0), 4),
        'Recall': round(recall_score(y_test, y_pred, zero_division=0), 4),
        'F1 Score': round(f1_score(y_test, y_pred, zero_division=0), 4),
        'NON_HEALTHY': int(np.bincount(y_test)[0]),
        'HEALTHY': int(np.bincount(y_test)[1] if len(np.bincount(y_test)) > 1 else 0),
        'Features': int(X_test.shape[1]),
        'Samples': int(len(y_test))
    }
    try:
        if y_prob is not None and len(np.unique(y_test)) > 1:
            metrics['AUC Score'] = round(roc_auc_score(y_test, y_prob, average='weighted'), 4)
            metrics['AUPR'] = round(average_precision_score(y_test, y_prob), 4)
        elif y_prob is None and len(np.unique(y_test)) > 1:
            metrics['AUC Score'] = round(roc_auc_score(y_test, y_pred, average='weighted'), 4)
            metrics['AUPR'] = None
        else:
            metrics['AUC Score'] = None
            metrics['AUPR'] = None
    except ValueError:
        metrics['AUC Score'] = None
        metrics['AUPR'] = None
    return metrics


def evaluate_feature_importance(models, X, y, feature_selection_method, dataset_type, output_dir, bioproject_id=None):
    """
    Compute per-model feature importance:
      - LogReg / SVM-Linear: model coefficients
      - SVM-RBF: permutation importance (F1-scored, 5 repeats)
      - RF: built-in impurity-based feature_importances_
    """
    feature_importance_list = []
    for model_name, model in models.items():
        print(f"Processing model: {model_name}, type: {type(model)}")
        importance_df = None

        if len(np.unique(y)) < 2:
            print(f"Skipping feature importance for {model_name} in {feature_selection_method} (Dataset: {dataset_type}): only one class present")
            continue

        if 'LogReg' in model_name:
            try:
                importance_df = pd.DataFrame({
                    'Feature': X.columns, 'Importance': model.coef_.flatten(), 'Model': model_name,
                    'Feature_Selection': feature_selection_method, 'Prevalence': prevalence, 'Dataset': dataset_type
                })
            except AttributeError as e:
                print(f"Error computing coef_ for {model_name}: {e}")
                continue

        elif 'SVM_RBF' in model_name:
            try:
                perm_importance = permutation_importance(model, X, y, scoring='f1', n_repeats=5, random_state=42, n_jobs=-1)
                importance_df = pd.DataFrame({
                    'Feature': X.columns, 'Importance': perm_importance.importances_mean, 'Model': model_name,
                    'Feature_Selection': feature_selection_method, 'Prevalence': prevalence, 'Dataset': dataset_type
                })
            except Exception as e:
                print(f"Error computing permutation_importance for {model_name}: {e}")
                continue

        elif 'SVM_LIN' in model_name:
            try:
                importance_df = pd.DataFrame({
                    'Feature': X.columns, 'Importance': model.coef_.flatten(), 'Model': model_name,
                    'Feature_Selection': feature_selection_method, 'Prevalence': prevalence, 'Dataset': dataset_type
                })
            except AttributeError as e:
                print(f"Error computing coef_ for {model_name}: {e}")
                continue

        elif 'RF' in model_name:
            try:
                importance_df = pd.DataFrame({
                    'Feature': X.columns, 'Importance': model.feature_importances_, 'Model': model_name,
                    'Feature_Selection': feature_selection_method, 'Prevalence': prevalence, 'Dataset': dataset_type
                })
            except AttributeError as e:
                print(f"Error computing feature_importances_ for {model_name}: {e}")
                continue
        else:
            print(f"Skipping unknown model type: {model_name}")
            continue

        if bioproject_id:
            importance_df['BioProject'] = bioproject_id

        if importance_df is not None:
            output_file = os.path.join(
                output_dir,
                f'feature_importance_{model_name}_{feature_selection_method}_'
                f'{"bioproject_" + bioproject_id + "_" if bioproject_id else ""}{dataset_type}_{prevalence}.tsv'
            )
            importance_df.to_csv(output_file, sep='\t', index=False)
            print(f"Feature importance saved to: {output_file}")
            feature_importance_list.append(importance_df)
        else:
            print(f"No feature importance generated for {model_name} in {feature_selection_method} (Dataset: {dataset_type})")

    return feature_importance_list


def plot_roc_curves(data_type, models_all_features, X_validation, y_validation, models_pfi, X_validation_pfi,
                     y_validation_pfi, models_rfe, X_validation_rfe, y_validation_rfe, models_lasso,
                     X_validation_lasso, y_validation_lasso):
    """2x2 ROC-curve panel (All Features / PFI / RFECV / LASSO) for one dataset (test or validation)."""
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_context('paper', font_scale=1.5)
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    axes = axes.flatten()

    panels = [
        (axes[0], models_all_features, X_validation, y_validation, 'All Features'),
        (axes[1], models_pfi, X_validation_pfi, y_validation_pfi, 'Permutation Feature Importance'),
        (axes[2], models_rfe, X_validation_rfe, y_validation_rfe, 'Recursive Feature Elimination'),
        (axes[3], models_lasso, X_validation_lasso, y_validation_lasso, 'LASSO-based Features'),
    ]

    for ax, models, X_set, y_set, panel_title in panels:
        for model_name, model in models.items():
            y_prob = model.predict_proba(X_set)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_set)
            fpr, tpr, _ = roc_curve(y_set, y_prob)
            roc_auc = auc(fpr, tpr)
            ax.plot(fpr, tpr, label=f'{model_name} AUC = {roc_auc:.2f}', lw=2)
        ax.plot([0, 1], [0, 1], linestyle='--', color='gray')
        ax.set_title(f'ROC Curve ({panel_title})', fontsize=14, fontweight='bold')
        ax.set_xlabel('False Positive Rate', fontsize=12, fontweight='bold')
        ax.set_ylabel('True Positive Rate', fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=10, prop={'weight': 'bold'})
        ax.grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.suptitle(f'ROC Curves for {data_type} Dataset (Prevalence: {prevalence})', fontsize=16, fontweight='bold')
    plt.savefig(os.path.join(output_dir, f'ROC_CURVE_{data_type}_dataset_{prevalence}.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_roc_curves_bioproject(models, X_test, y_test, project_id, phenotype_counts_str, output_dir):
    """Compute and save ROC curve data for all models, restricted to a single BioProject."""
    os.makedirs(output_dir, exist_ok=True)

    if len(np.unique(y_test)) < 2:
        print(f"Skipping ROC curve for BioProject {project_id}: only one class present in y_test")
        return pd.DataFrame()

    roc_data_list = []
    for model_name, model in models.items():
        try:
            y_prob = model.predict_proba(X_test)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_test)
            fpr, tpr, thresholds = roc_curve(y_test, y_prob)
            roc_auc = auc(fpr, tpr)
            roc_data = pd.DataFrame({
                'FPR': fpr, 'TPR': tpr, 'Thresholds': thresholds, 'Model': model_name,
                'BioProject': project_id, 'Phenotype_Counts': phenotype_counts_str, 'ROC_AUC': roc_auc
            })
            roc_data_list.append(roc_data)
            roc_data.to_csv(os.path.join(output_dir, f'roc_data_{project_id}_{model_name}_{prevalence}.tsv'), sep='\t', index=False)
        except ValueError as e:
            print(f"Error computing ROC for {project_id}, {model_name}: {e}")
            continue

    return pd.concat(roc_data_list, ignore_index=True) if roc_data_list else pd.DataFrame()


def save_roc_data(models, X_validation, y_validation, feature_selection_method):
    """Collect raw ROC curve points (FPR/TPR/thresholds) for all models under one feature-selection method."""
    roc_data = []
    for model_name, model in models.items():
        if len(np.unique(y_validation)) < 2:
            continue
        y_prob = model.predict_proba(X_validation)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_validation)
        fpr, tpr, thresholds = roc_curve(y_validation, y_prob)
        roc_auc = auc(fpr, tpr)
        for fp, tp, thresh in zip(fpr, tpr, thresholds):
            roc_data.append({
                "Model": model_name, "Feature_Selection": feature_selection_method,
                "FPR": fp, "TPR": tp, "Thresholds": thresh, "AUC": roc_auc
            })
    return roc_data


def plot_precision_recall_curves(data_type, models_all_features, X_validation, y_validation, models_pfi,
                                  X_validation_pfi, y_validation_pfi, models_rfe, X_validation_rfe,
                                  y_validation_rfe, models_lasso, X_validation_lasso, y_validation_lasso):
    """2x2 Precision-Recall curve panel (All Features / PFI / RFECV / LASSO) for one dataset."""
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_context('paper', font_scale=1.5)
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    axes = axes.flatten()

    panels = [
        (axes[0], models_all_features, X_validation, y_validation, 'All Features'),
        (axes[1], models_pfi, X_validation_pfi, y_validation_pfi, 'PFI'),
        (axes[2], models_rfe, X_validation_rfe, y_validation_rfe, 'RFE'),
        (axes[3], models_lasso, X_validation_lasso, y_validation_lasso, 'Lasso'),
    ]

    for ax, models, X_set, y_set, panel_title in panels:
        for model_name, model in models.items():
            y_prob = model.predict_proba(X_set)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_set)
            precision, recall, _ = precision_recall_curve(y_set, y_prob)
            ap = average_precision_score(y_set, y_prob)
            ax.plot(recall, precision, label=f'{model_name} AP = {ap:.2f}', lw=2)
        ax.set_title(f'Precision-Recall Curve ({panel_title})', fontsize=14)
        ax.set_xlabel('Recall', fontsize=12)
        ax.set_ylabel('Precision', fontsize=12)
        ax.legend(loc='best', fontsize=10)
        ax.grid(True)

    plt.tight_layout()
    plt.suptitle(f'Precision-Recall Curves for {data_type} dataset ({prevalence})', fontsize=15, y=1.05)
    plt.savefig(os.path.join(output_dir, f'PR_CURVE_{data_type}_dataset_{prevalence}.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)


def save_precision_recall_data(models, X_validation, y_validation, feature_selection_method):
    """Collect raw Precision-Recall curve points for all models under one feature-selection method."""
    pr_data = []
    for model_name, model in models.items():
        if len(np.unique(y_validation)) < 2:
            continue
        y_prob = model.predict_proba(X_validation)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_validation)
        precision, recall, _ = precision_recall_curve(y_validation, y_prob)
        ap = average_precision_score(y_validation, y_prob)
        for p, r in zip(precision, recall):
            pr_data.append({
                "Model": model_name, "Feature_Selection": feature_selection_method,
                "Precision": p, "Recall": r, "Average_Precision": ap
            })
    return pr_data


def compute_confusion_matrix(models, X, y, sample_ids, feature_selection_method, dataset_type, output_dir,
                              prevalence=prevalence, bioproject_id=None):
    """Compute and save the raw confusion matrix (TN/FP/FN/TP counts) for each model."""
    cm_results = []
    cm_dir = os.path.join(output_dir, 'confusion_matrix_results')
    os.makedirs(cm_dir, exist_ok=True)

    for model_name, model in models.items():
        if len(X) != len(y):
            print(f"Warning: Inconsistent sample counts for {model_name} ({feature_selection_method}, {dataset_type}): X has {len(X)} samples, y has {len(y)} samples")
            continue
        try:
            y_pred = model.predict(X)
            if len(y_pred) != len(y):
                print(f"Warning: Prediction length mismatch for {model_name} ({feature_selection_method}, {dataset_type}): y_pred has {len(y_pred)} samples, y has {len(y)} samples")
                continue

            cm = confusion_matrix(y, y_pred, labels=[0, 1])
            cm_df = pd.DataFrame({
                'Metric': ['True_Negative', 'False_Positive', 'False_Negative', 'True_Positive'],
                'Count': [cm[0, 0], cm[0, 1], cm[1, 0], cm[1, 1]],
                'Model': model_name, 'Feature_Selection': feature_selection_method,
                'Prevalence': prevalence, 'Dataset': dataset_type
            })
            if bioproject_id:
                cm_df['BioProject'] = bioproject_id
                output_file = os.path.join(cm_dir, f'confusion_matrix_{model_name}_{feature_selection_method}_{bioproject_id}_{prevalence}.tsv')
            else:
                output_file = os.path.join(cm_dir, f'confusion_matrix_{model_name}_{feature_selection_method}_{dataset_type}_{prevalence}.tsv')

            cm_df.to_csv(output_file, sep='\t', index=False)
            print(f"Confusion matrix saved to: {output_file}")
            cm_results.append(cm_df)
        except ValueError as e:
            print(f"Error computing confusion matrix for {model_name} ({feature_selection_method}, {dataset_type}): {e}")
            continue

    return cm_results


def sample_classification_details(models, X, y, sample_ids, feature_selection_method, dataset_type, output_dir,
                                    prevalence=prevalence, bioproject_id=None):
    """Save per-sample true/predicted labels and correctness status for each model."""
    classification_results = []
    cm_dir = os.path.join(output_dir, 'confusion_matrix_results')
    os.makedirs(cm_dir, exist_ok=True)

    for model_name, model in models.items():
        y_pred = model.predict(X)
        classification_status = (y_pred == y).astype(int)
        classification_df = pd.DataFrame({
            'Run_ID': sample_ids.get('Run_ID', pd.Series([None] * len(y))),
            'BioProject': sample_ids.get('BioProject', pd.Series([None] * len(y))),
            'Phenotype': sample_ids.get('Phenotype', pd.Series([None] * len(y))),
            'True_Label': y, 'Predicted_Label': y_pred,
            'Classification_Status': ['Correct' if s == 1 else 'Incorrect' for s in classification_status],
            'Model': model_name, 'Feature_Selection': feature_selection_method,
            'Prevalence': prevalence, 'Dataset': dataset_type
        })
        if bioproject_id:
            classification_df['BioProject'] = bioproject_id
            output_file = os.path.join(cm_dir, f'classification_details_{model_name}_{feature_selection_method}_bioproject_{bioproject_id}_{prevalence}.tsv')
        else:
            output_file = os.path.join(cm_dir, f'classification_details_{model_name}_{feature_selection_method}_{dataset_type}_{prevalence}.tsv')

        classification_df.to_csv(output_file, sep='\t', index=False)
        print(f"Classification details saved to: {output_file}")
        classification_results.append(classification_df)

    return classification_results


def plot_confusion_matrix(models, X, y, sample_ids, feature_selection_method, dataset_type, output_dir,
                           prevalence, bioproject_id=None):
    """Save the confusion matrix table (counts + percentages) and render an annotated heatmap per model."""
    cm_results = []
    cm_dir = os.path.join(output_dir, 'confusion_matrix_results')
    os.makedirs(cm_dir, exist_ok=True)

    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_context('paper', font_scale=1.8)

    for model_name, model in models.items():
        y_pred = model.predict(X)
        cm = confusion_matrix(y, y_pred, labels=[0, 1])
        total_samples = len(y)
        cm_percent = (cm / total_samples * 100).round(2)

        cm_df = pd.DataFrame({
            'Metric': ['True_Negative', 'False_Positive', 'False_Negative', 'True_Positive'],
            'Count': [cm[0, 0], cm[0, 1], cm[1, 0], cm[1, 1]],
            'Percentage': [cm_percent[0, 0], cm_percent[0, 1], cm_percent[1, 0], cm_percent[1, 1]],
            'Model': model_name, 'Feature_Selection': feature_selection_method,
            'Prevalence': prevalence, 'Dataset': dataset_type
        })
        if bioproject_id:
            cm_df['BioProject'] = bioproject_id
        cm_results.append(cm_df)

        fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
        sns.heatmap(cm, annot=False, cmap='Blues', cbar=True,
                    xticklabels=['Non-Healthy (0)', 'Healthy (1)'],
                    yticklabels=['Non-Healthy (0)', 'Healthy (1)'], ax=ax)
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                text = f'{int(cm[i, j])}\n({cm_percent[i, j]:.2f}%)'
                ax.text(j + 0.5, i + 0.5, text, ha='center', va='center',
                        fontsize=14, color='black', weight='bold',
                        bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white', alpha=0.7))
        ax.set_xlabel('Predicted Label', fontsize=14, weight='bold')
        ax.set_ylabel('True Label', fontsize=14, weight='bold')

        title = f'Confusion Matrix: {model_name} ({feature_selection_method})\nDataset: {dataset_type}, Prevalence: {prevalence}'
        if bioproject_id:
            title += f'\nBioProject: {bioproject_id}'
            phenotype_counts = sample_ids['Phenotype'].value_counts().to_dict()
            title += f'\nPhenotypes: {", ".join([f"{k} ({v})" for k, v in phenotype_counts.items()])}'
        title += f'\nSamples: {total_samples}'
        ax.set_title(title, fontsize=16, pad=20, weight='bold')

        cbar = ax.collections[0].colorbar
        cbar.set_label('Count', fontsize=12, weight='bold')
        cbar.ax.tick_params(labelsize=10)
        plt.tight_layout()

        if bioproject_id:
            output_file = os.path.join(cm_dir, f'confusion_matrix_plot_{model_name}_{feature_selection_method}_bioproject_{bioproject_id}_{prevalence}.png')
        else:
            output_file = os.path.join(cm_dir, f'confusion_matrix_plot_{model_name}_{feature_selection_method}_{dataset_type}_{prevalence}.png')

        plt.savefig(output_file, dpi=600, bbox_inches='tight')
        plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
        print(f"Confusion matrix plot saved to: {output_file} and {output_file.replace('.png', '.pdf')}")
        plt.close(fig)

    return cm_results


def analyze_misclassifications(models, X, y, sample_ids, feature_selection_method, dataset_type, output_dir, prevalence):
    """Identify, save, and plot misclassified samples by phenotype, per model."""
    cm_dir = os.path.join(output_dir, 'confusion_matrix_results')
    os.makedirs(cm_dir, exist_ok=True)
    misclass_df_list = []

    for model_name, model in models.items():
        try:
            y_pred = model.predict(X)
            if len(y_pred) != len(y):
                print(f"Warning: Prediction length mismatch for {model_name} ({feature_selection_method}, {dataset_type}): y_pred has {len(y_pred)} samples, y has {len(y)} samples")
                continue

            misclassified = y_pred != y
            if not misclassified.any():
                print(f"No misclassifications for {model_name} ({feature_selection_method}, {dataset_type})")
                continue

            misclass_df = sample_ids.loc[misclassified].copy()
            misclass_df['Model'] = model_name
            misclass_df['Feature_Selection'] = feature_selection_method
            misclass_df['True_Label'] = y[misclassified]
            misclass_df['Predicted_Label'] = y_pred[misclassified]
            misclass_df['Prevalence'] = prevalence
            misclass_df['Dataset'] = dataset_type

            output_file = os.path.join(cm_dir, f'misclassifications_{model_name}_{feature_selection_method}_{dataset_type}_{prevalence}.tsv')
            misclass_df.to_csv(output_file, sep='\t', index=False)
            print(f"Misclassification data saved to: {output_file}")
            misclass_df_list.append(misclass_df)

            plt.style.use('seaborn-v0_8-whitegrid')
            sns.set_context('paper', font_scale=1.2)
            plt.figure(figsize=(10, 6))
            sns.countplot(data=misclass_df, x='Phenotype', hue='Model', palette='Set2')
            plt.title(f'Misclassifications by Phenotype: {model_name} ({feature_selection_method}, {dataset_type}, Prevalence: {prevalence})', fontsize=12, pad=15)
            plt.xlabel('Phenotype', fontsize=10)
            plt.ylabel('Count', fontsize=10)
            plt.xticks(rotation=45, ha='right')
            plt.legend(title='Model', fontsize=8)
            plt.tight_layout()

            plot_file = os.path.join(cm_dir, f'misclass_plot_{model_name}_{feature_selection_method}_{dataset_type}_{prevalence}.png')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            print(f"Misclassification plot saved to: {plot_file}")
            plt.close()

        except Exception as e:
            print(f"Error analyzing misclassifications for {model_name} ({feature_selection_method}, {dataset_type}): {e}")
            continue

    return misclass_df_list


def generate_classification_report(models, X, y, sample_ids, feature_selection_method, dataset_type, output_dir,
                                     prevalence, bioproject_id=None):
    """Compute and save the sklearn per-class classification report (precision/recall/F1) for each model."""
    report_results = []
    report_dir = os.path.join(output_dir, 'classification_reports')
    os.makedirs(report_dir, exist_ok=True)

    for model_name, model in models.items():
        try:
            if len(X) != len(y):
                print(f"Warning: Inconsistent sample counts for {model_name} ({feature_selection_method}, {dataset_type}): X has {len(X)} samples, y has {len(y)} samples")
                continue

            y_pred = model.predict(X)
            report_dict = classification_report(y, y_pred, output_dict=True, zero_division=0)

            report_data = []
            for label, metrics in report_dict.items():
                if isinstance(metrics, dict):
                    report_data.append({
                        'Model': model_name, 'Feature_Selection': feature_selection_method, 'Class_Label': label,
                        'Precision': round(metrics['precision'], 4), 'Recall': round(metrics['recall'], 4),
                        'F1_Score': round(metrics['f1-score'], 4), 'Support': int(metrics['support']),
                        'Prevalence': prevalence, 'Dataset': dataset_type
                    })

            report_df = pd.DataFrame(report_data)
            if bioproject_id:
                report_df['BioProject'] = bioproject_id
                output_file = os.path.join(report_dir, f'classification_report_{model_name}_{feature_selection_method}_bioproject_{bioproject_id}_{prevalence}.tsv')
            else:
                output_file = os.path.join(report_dir, f'classification_report_{model_name}_{feature_selection_method}_{dataset_type}_{prevalence}.tsv')

            report_df.to_csv(output_file, sep='\t', index=False)
            print(f"Classification report saved to: {output_file}")
            report_results.append(report_df)

        except Exception as e:
            print(f"Error generating classification report for {model_name} ({feature_selection_method}, {dataset_type}): {e}")
            continue

    return report_results


def plot_classification_report(models, X, y, sample_ids, feature_selection_method, dataset_type, output_dir,
                                 prevalence, bioproject_id=None):
    """Render a grouped bar chart (Precision/Recall/F1 per class) for each model's classification report."""
    report_results = []
    report_dir = os.path.join(output_dir, 'classification_reports')
    os.makedirs(report_dir, exist_ok=True)

    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_context('paper', font_scale=1.8)

    for model_name, model in models.items():
        try:
            y_pred = model.predict(X)
            report_dict = classification_report(y, y_pred, output_dict=True, zero_division=0)

            report_data = []
            for label, metrics in report_dict.items():
                if isinstance(metrics, dict):
                    report_data.append({
                        'Class_Label': label, 'Precision': metrics['precision'], 'Recall': metrics['recall'],
                        'F1_Score': metrics['f1-score'], 'Support': metrics['support']
                    })
            report_df = pd.DataFrame(report_data)
            report_df['Model'] = model_name
            report_df['Feature_Selection'] = feature_selection_method
            report_df['Prevalence'] = prevalence
            report_df['Dataset'] = dataset_type
            if bioproject_id:
                report_df['BioProject'] = bioproject_id
            report_results.append(report_df)

            plot_df = report_df.melt(id_vars=['Class_Label'], value_vars=['Precision', 'Recall', 'F1_Score'],
                                      var_name='Metric', value_name='Score')
            plt.figure(figsize=(12, 8), dpi=300)
            sns.barplot(data=plot_df, x='Class_Label', y='Score', hue='Metric', palette='Set2')

            phenotype_counts = sample_ids['Phenotype'].value_counts().to_dict()
            phenotype_title = ", ".join([f"{key} ({value})" for key, value in phenotype_counts.items()])
            title = f'Classification Report: {model_name} ({feature_selection_method})\nDataset: {dataset_type}, Prevalence: {prevalence}'
            if bioproject_id:
                title += f'\nBioProject: {bioproject_id}\nPhenotypes: {phenotype_title}'
            plt.title(title, fontsize=16, pad=20, weight='bold')
            plt.xlabel('Class Label', fontsize=14, weight='bold')
            plt.ylabel('Score', fontsize=14, weight='bold')
            plt.legend(title='Metric', fontsize=12, frameon=True, edgecolor='black')
            plt.ylim(0, 1.05)
            plt.tight_layout()

            if bioproject_id:
                plot_file = os.path.join(report_dir, f'classification_report_plot_{model_name}_{feature_selection_method}_{bioproject_id}_{prevalence}.png')
            else:
                plot_file = os.path.join(report_dir, f'classification_report_plot_{model_name}_{feature_selection_method}_{dataset_type}_{prevalence}.png')

            plt.savefig(plot_file, dpi=600, bbox_inches='tight')
            plt.savefig(plot_file.replace('.png', '.pdf'), bbox_inches='tight')
            print(f"Classification report plot saved to: {plot_file} and {plot_file.replace('.png', '.pdf')}")
            plt.close()

        except Exception as e:
            print(f"Error plotting classification report for {model_name} ({feature_selection_method}, {dataset_type}): {e}")
            continue

    return report_results


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

sample_ids_validation = copy_df[['Run_ID', 'BioProject', 'Phenotype']].loc[X_validation.index]
sample_ids_test = df_copy[['Run_ID', 'BioProject', 'Phenotype']].loc[X_test.index]

print("\nEvaluating feature importance on test dataset...")
X_test_rfe = X_test[[col for col in rfe_features if col in X_test.columns]]
X_test_pfi = X_test[[col for col in pfi_features if col in X_test.columns]]
X_test_lasso = X_test[[col for col in lasso_features if col in X_test.columns]]


# ------------------------------------------------------------------------------
# Step 1: Evaluate all models on the test dataset (per feature-selection method)
# ------------------------------------------------------------------------------
print("\nEvaluating models on test dataset...")
results_test = []
for model_name, model in models_all_features.items():
    matrix = evaluate_model(model, X_test, y_test, model.predict_proba(X_test)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_test))
    results_test.append({'Feature_Selection': 'All_Features', 'Model': model_name, **matrix})
    print(f"{model_name}\t {matrix}")

for model_name, model in models_rfe.items():
    matrix = evaluate_model(model, X_test_rfe, y_test, model.predict_proba(X_test_rfe)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_test_rfe))
    results_test.append({'Feature_Selection': 'RFE', 'Model': model_name, **matrix})
    print(f"{model_name}\t {matrix}")

for model_name, model in models_pfi.items():
    matrix = evaluate_model(model, X_test_pfi, y_test, model.predict_proba(X_test_pfi)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_test_pfi))
    results_test.append({'Feature_Selection': 'PFI', 'Model': model_name, **matrix})
    print(f"{model_name}\t {matrix}")

for model_name, model in models_lasso.items():
    matrix = evaluate_model(model, X_test_lasso, y_test, model.predict_proba(X_test_lasso)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_test_lasso))
    results_test.append({'Feature_Selection': 'LASSO', 'Model': model_name, **matrix})
    print(f"{model_name}\t {matrix}")

pd.DataFrame(results_test).to_csv(os.path.join(output_dir, f'Performance_on_test_dataset_{prevalence}.tsv'), sep='\t', index=False)


# ------------------------------------------------------------------------------
# Step 2: Evaluate all models on the independent validation dataset
# ------------------------------------------------------------------------------
print("\nEvaluating models on validation dataset...")
results_validation = []
for model_name, model in models_all_features.items():
    matrix = evaluate_model(model, X_validation, y_validation, model.predict_proba(X_validation)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_validation))
    results_validation.append({'Feature_Selection': 'All_Features', 'Model': model_name, **matrix})
    print(f"{model_name}\t {matrix}")

for model_name, model in models_rfe.items():
    matrix = evaluate_model(model, X_validation_rfe, y_validation_rfe, model.predict_proba(X_validation_rfe)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_validation_rfe))
    results_validation.append({'Feature_Selection': 'RFE', 'Model': model_name, **matrix})
    print(f"{model_name}\t {matrix}")

for model_name, model in models_pfi.items():
    matrix = evaluate_model(model, X_validation_pfi, y_validation_pfi, model.predict_proba(X_validation_pfi)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_validation_pfi))
    results_validation.append({'Feature_Selection': 'PFI', 'Model': model_name, **matrix})
    print(f"{model_name}\t {matrix}")

for model_name, model in models_lasso.items():
    matrix = evaluate_model(model, X_validation_lasso, y_validation_lasso, model.predict_proba(X_validation_lasso)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_validation_lasso))
    results_validation.append({'Feature_Selection': 'LASSO', 'Model': model_name, **matrix})
    print(f"{model_name}\t {matrix}")

pd.DataFrame(results_validation).to_csv(os.path.join(output_dir, f'Performance_on_validation_dataset_{prevalence}.tsv'), sep='\t', index=False)


# ------------------------------------------------------------------------------
# Step 3: Per-class classification reports (test dataset)
# ------------------------------------------------------------------------------
print("\nGenerating classification reports for test dataset...")
all_classification_reports_test = []
for models, X_set, method in [(models_all_features, X_test, "All_Features"), (models_rfe, X_test_rfe, "RFE"),
                               (models_pfi, X_test_pfi, "PFI"), (models_lasso, X_test_lasso, "LASSO")]:
    all_classification_reports_test.extend(generate_classification_report(models, X_set, y_test, sample_ids_test, method, "test", output_dir, prevalence))
    all_classification_reports_test.extend(plot_classification_report(models, X_set, y_test, sample_ids_test, method, "test", output_dir, prevalence))

if all_classification_reports_test:
    out_path = os.path.join(output_dir, 'classification_reports', f'classification_reports_test_dataset_{prevalence}.tsv')
    pd.concat(all_classification_reports_test, ignore_index=True).to_csv(out_path, sep='\t', index=False)
    print(f"Combined classification report results saved to: {out_path}")


# ------------------------------------------------------------------------------
# Step 4: Per-class classification reports (validation dataset)
# ------------------------------------------------------------------------------
print("\nGenerating classification reports for validation dataset...")
all_classification_reports_validation = []
for models, X_set, y_set, method in [
    (models_all_features, X_validation, y_validation, "All_Features"),
    (models_rfe, X_validation_rfe, y_validation_rfe, "RFE"),
    (models_pfi, X_validation_pfi, y_validation_pfi, "PFI"),
    (models_lasso, X_validation_lasso, y_validation_lasso, "LASSO")
]:
    all_classification_reports_validation.extend(generate_classification_report(models, X_set, y_set, sample_ids_validation, method, "validation", output_dir, prevalence))
    all_classification_reports_validation.extend(plot_classification_report(models, X_set, y_set, sample_ids_validation, method, "validation", output_dir, prevalence))

if all_classification_reports_validation:
    out_path = os.path.join(output_dir, 'classification_reports', f'classification_reports_validation_dataset_{prevalence}.tsv')
    pd.concat(all_classification_reports_validation, ignore_index=True).to_csv(out_path, sep='\t', index=False)
    print(f"Combined classification report results saved to: {out_path}")


# ------------------------------------------------------------------------------
# Step 5: ROC and Precision-Recall curves (test and validation datasets)
# ------------------------------------------------------------------------------
print("\nPlotting ROC and Precision-Recall curves for test datasets...")
plot_roc_curves("test", models_all_features, X_test, y_test, models_pfi, X_test_pfi, y_test, models_rfe, X_test_rfe, y_test, models_lasso, X_test_lasso, y_test)
plot_precision_recall_curves("test", models_all_features, X_test, y_test, models_pfi, X_test_pfi, y_test, models_rfe, X_test_rfe, y_test, models_lasso, X_test_lasso, y_test)

print("\nPlotting ROC and Precision-Recall curves for validation datasets...")
plot_roc_curves("validation", models_all_features, X_validation, y_validation, models_pfi, X_validation_pfi, y_validation_pfi, models_rfe, X_validation_rfe, y_validation_rfe, models_lasso, X_validation_lasso, y_validation_lasso)
plot_precision_recall_curves("validation", models_all_features, X_validation, y_validation, models_pfi, X_validation_pfi, y_validation_pfi, models_rfe, X_validation_rfe, y_validation_rfe, models_lasso, X_validation_lasso, y_validation_lasso)


# ------------------------------------------------------------------------------
# Step 6: Save raw ROC / Precision-Recall curve data (test dataset)
# ------------------------------------------------------------------------------
all_roc_data = []
all_roc_data.extend(save_roc_data(models_all_features, X_test, y_test, "All_Features"))
all_roc_data.extend(save_roc_data(models_pfi, X_test_pfi, y_test, "PFI"))
all_roc_data.extend(save_roc_data(models_rfe, X_test_rfe, y_test, "RFE"))
all_roc_data.extend(save_roc_data(models_lasso, X_test_lasso, y_test, "LASSO"))
roc_df = pd.DataFrame(all_roc_data)
roc_out_path = os.path.join(roc_data_dir, f'roc_data_test_dataset_{prevalence}.csv')
roc_df.to_csv(roc_out_path, index=False)
print(f"ROC data saved to: {roc_out_path}")

all_pr_data = []
all_pr_data.extend(save_precision_recall_data(models_all_features, X_test, y_test, "All_Features"))
all_pr_data.extend(save_precision_recall_data(models_pfi, X_test_pfi, y_test, "PFI"))
all_pr_data.extend(save_precision_recall_data(models_rfe, X_test_rfe, y_test, "RFE"))
all_pr_data.extend(save_precision_recall_data(models_lasso, X_test_lasso, y_test, "LASSO"))
pr_df_test = pd.DataFrame(all_pr_data)
pr_out_path = os.path.join(roc_data_dir, f'precision_recall_for_test_dataset_{prevalence}.tsv')
pr_df_test.to_csv(pr_out_path, sep='\t', index=False)
print(f"Precision-Recall data saved to: {pr_out_path}")


# ------------------------------------------------------------------------------
# Step 7: Confusion matrices, classification details, misclassifications (test)
# ------------------------------------------------------------------------------
print("\nEvaluating confusion matrix, classification details, and misclassifications on test dataset...")
sample_ids_test = df_copy[['Run_ID', 'BioProject', 'Phenotype']].loc[X_test.index]
all_cm_results, all_classification_results, all_misclass_results = [], [], []

for models, X_set, method in [(models_all_features, X_test, "All_Features"), (models_rfe, X_test_rfe, "RFE"),
                               (models_pfi, X_test_pfi, "PFI"), (models_lasso, X_test_lasso, "LASSO")]:
    all_cm_results.extend(compute_confusion_matrix(models, X_set, y_test, sample_ids_test, method, "test", output_dir, prevalence))
    all_classification_results.extend(sample_classification_details(models, X_set, y_test, sample_ids_test, method, "test", output_dir, prevalence))
    all_cm_results.extend(plot_confusion_matrix(models, X_set, y_test, sample_ids_test, method, "test", output_dir, prevalence))
    all_misclass_results.extend(analyze_misclassifications(models, X_set, y_test, sample_ids_test, method, "test", output_dir, prevalence))

pd.concat(all_cm_results, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'confusion_matrix_test_dataset_{prevalence}.tsv'), sep='\t', index=False)
pd.concat(all_classification_results, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'classification_details_test_dataset_{prevalence}.tsv'), sep='\t', index=False)
if all_misclass_results:
    pd.concat(all_misclass_results, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'misclassifications_test_dataset_{prevalence}.tsv'), sep='\t', index=False)
print("Test dataset confusion matrix / classification detail / misclassification results saved.")


# ------------------------------------------------------------------------------
# Step 8: Confusion matrices, classification details, misclassifications (validation)
# ------------------------------------------------------------------------------
print("\nEvaluating confusion matrix, classification details, and misclassifications on validation dataset...")
sample_ids_validation = copy_df[['Run_ID', 'BioProject', 'Phenotype']].loc[X_validation.index]
all_cm_results, all_classification_results, all_misclass_results = [], [], []

for models, X_set, y_set, method in [
    (models_all_features, X_validation, y_validation, "All_Features"),
    (models_rfe, X_validation_rfe, y_validation_rfe, "RFE"),
    (models_pfi, X_validation_pfi, y_validation_pfi, "PFI"),
    (models_lasso, X_validation_lasso, y_validation_lasso, "LASSO")
]:
    all_cm_results.extend(compute_confusion_matrix(models, X_set, y_set, sample_ids_validation, method, "validation", output_dir, prevalence))
    all_classification_results.extend(sample_classification_details(models, X_set, y_set, sample_ids_validation, method, "validation", output_dir, prevalence))
    all_cm_results.extend(plot_confusion_matrix(models, X_set, y_set, sample_ids_validation, method, "validation", output_dir, prevalence))
    all_misclass_results.extend(analyze_misclassifications(models, X_set, y_set, sample_ids_validation, method, "validation", output_dir, prevalence))

pd.concat(all_cm_results, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'confusion_matrix_validation_dataset_{prevalence}.tsv'), sep='\t', index=False)
pd.concat(all_classification_results, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'classification_details_validation_dataset_{prevalence}.tsv'), sep='\t', index=False)
if all_misclass_results:
    pd.concat(all_misclass_results, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'misclassifications_validation_dataset_{prevalence}.tsv'), sep='\t', index=False)
print("Validation dataset confusion matrix / classification detail / misclassification results saved.")


# ------------------------------------------------------------------------------
# Step 9: Per-BioProject evaluation - test dataset
# ------------------------------------------------------------------------------
print("\nEvaluating models on each BioProject in test datasets...")
test_grouped_data = df_copy.loc[test_indices].groupby('BioProject')

all_metrics_test = []
all_bioproject_feature_importance_test = []
all_cm_results_test, all_classification_results_test, all_misclass_results_test = [], [], []
all_roc_data_test = []
bioproject_to_phenotype = {}

for project_id, group in test_grouped_data:
    phenotype = group["Phenotype"].unique()
    phenotype_str = ",".join(phenotype)
    phenotype_counts = group["Phenotype"].value_counts().to_dict()
    phenotype_counts_str = f"{project_id}: {','.join([f'{key}({value})' for key, value in phenotype_counts.items()])}"
    bioproject_to_phenotype[project_id] = phenotype_counts_str

    test_df_bioproject = group.drop(columns=[col for col in columns_to_drop if col in group.columns])
    y_test_bioproject = test_df_bioproject['Class_Label'].astype(int)
    X_test_all = test_df_bioproject.drop(columns=['Class_Label'])

    X_test_bioproject = X_test_all[[col for col in X_train.columns if col in X_test_all.columns]]
    X_test_rfe = X_test_all[[col for col in rfe_features if col in X_test_all.columns]].dropna()
    X_test_pfi = X_test_all[[col for col in pfi_features if col in X_test_all.columns]]
    X_test_lasso = X_test_all[[col for col in lasso_features if col in X_test_all.columns]]

    y_test_rfe = y_test_bioproject.loc[X_test_rfe.index]
    y_test_pfi = y_test_bioproject.loc[X_test_pfi.index]
    y_test_lasso = y_test_bioproject.loc[X_test_lasso.index]
    sample_ids_bioproject = group[['Run_ID', 'BioProject', 'Phenotype']]

    print(f"BioProject {project_id} - y_test_bioproject unique values: {np.unique(y_test_bioproject)}")
    print(f"BioProject {project_id} - X_test_bioproject shape: {X_test_bioproject.shape}")
    print(f"BioProject {project_id} - X_test_rfe shape: {X_test_rfe.shape}")
    print(f"BioProject {project_id} - y_test_rfe shape: {y_test_rfe.shape}")

    for models, X_set, y_set, method in [
        (models_all_features, X_test_bioproject, y_test_bioproject, "All_Features"),
        (models_pfi, X_test_pfi, y_test_pfi, "PFI"),
        (models_rfe, X_test_rfe, y_test_rfe, "RFE"),
        (models_lasso, X_test_lasso, y_test_lasso, "LASSO")
    ]:
        for model_name, model in models.items():
            matrix = evaluate_model(model, X_set, y_set, model.predict_proba(X_set)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_set))
            matrix.update({'BioProject': project_id, 'Phenotype': phenotype_str, 'Model': model_name,
                            'Feature_Selection': method, 'BioProject_Phenotype_Counts': phenotype_counts_str})
            all_metrics_test.append(matrix)

        all_bioproject_feature_importance_test.extend(evaluate_feature_importance(models, X_set, y_set, method, "bioproject", featrue_importance_dir, project_id))
        all_cm_results_test.extend(compute_confusion_matrix(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence, project_id))
        all_classification_results_test.extend(sample_classification_details(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence, project_id))
        all_cm_results_test.extend(plot_confusion_matrix(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence, project_id))
        all_misclass_results_test.extend(analyze_misclassifications(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence))
        all_classification_reports_test.extend(generate_classification_report(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence, project_id))
        all_classification_reports_test.extend(plot_classification_report(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence, project_id))

        roc_data = plot_roc_curves_bioproject(models, X_set, y_set, project_id, phenotype_counts_str, roc_data_dir)
        if not roc_data.empty:
            roc_data['Feature_Selection'] = method
            all_roc_data_test.append(roc_data)

metrics_df_test = pd.DataFrame(all_metrics_test)
metrics_df_test.to_csv(os.path.join(output_dir, f'model_evaluation_on_each_bioproject_test_datasets_{prevalence}.csv'), index=False)
print(f"Test BioProject evaluation results saved to: {os.path.join(output_dir, f'model_evaluation_on_each_bioproject_test_datasets_{prevalence}.csv')}")

if all_bioproject_feature_importance_test:
    pd.concat(all_bioproject_feature_importance_test, ignore_index=True).to_csv(os.path.join(output_dir, f'feature_importance_bioproject_test_{prevalence}.tsv'), sep='\t', index=False)
if all_cm_results_test:
    pd.concat(all_cm_results_test, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'confusion_matrix_bioproject_test_{prevalence}.tsv'), sep='\t', index=False)
if all_classification_results_test:
    pd.concat(all_classification_results_test, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'classification_details_bioproject_test_{prevalence}.tsv'), sep='\t', index=False)
if all_misclass_results_test:
    pd.concat(all_misclass_results_test, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'misclassifications_bioproject_test_{prevalence}.tsv'), sep='\t', index=False)
if all_classification_reports_test:
    pd.concat(all_classification_reports_test, ignore_index=True).to_csv(os.path.join(output_dir, 'classification_reports', f'classification_reports_bioproject_test_{prevalence}.tsv'), sep='\t', index=False)

if all_roc_data_test:
    all_roc_data_test_df = pd.concat(all_roc_data_test, ignore_index=True)
    roc_test_path = os.path.join(roc_data_dir, f'roc_data_bioproject_test_{prevalence}.tsv')
    all_roc_data_test_df.to_csv(roc_test_path, sep='\t', index=False)
    print(f"Test ROC data saved to: {roc_test_path}")

    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_context('paper', font_scale=1.8)
    g = sns.FacetGrid(all_roc_data_test_df, col="BioProject", hue="Model", col_wrap=4, height=6, aspect=1.5)
    g.map(sns.lineplot, "FPR", "TPR", marker="o", linewidth=2.5, markersize=8)
    g.add_legend(title='Model', fontsize=12, frameon=True, edgecolor='black')
    g.set_axis_labels("False Positive Rate", "True Positive Rate", fontsize=14, weight='bold')
    for ax, project_id in zip(g.axes.flat, g.col_names):
        ax.set_title(bioproject_to_phenotype.get(project_id, project_id), fontsize=14, pad=10, weight='bold')
        ax.plot([0, 1], [0, 1], color='black', linestyle='--', alpha=0.7)
        ax.grid(True, linestyle='--', alpha=0.7)
    g.fig.tight_layout(pad=1.5)
    output_path = os.path.join(roc_data_dir, f'roc_facet_grid_test_{prevalence}.png')
    g.savefig(output_path, dpi=600, bbox_inches='tight')
    g.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"Test ROC facet grid saved to: {output_path} and {output_path.replace('.png', '.pdf')}")
    plt.close(g.fig)
else:
    print("No ROC data generated due to single-class BioProjects in test dataset")


# ------------------------------------------------------------------------------
# Step 10: Per-BioProject evaluation - validation dataset
# ------------------------------------------------------------------------------
print("\nEvaluating models on each BioProject in validation datasets...")
grouped_data = dataset.groupby('BioProject')

all_metrics = []
all_bioproject_feature_importance = []
all_cm_results, all_classification_results, all_misclass_results = [], [], []
all_roc_data = []
all_classification_reports = []
bioproject_to_phenotype = {}

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
    X_validation_rfe = X_validation_all[[col for col in rfe_features if col in X_validation_all.columns]].dropna()
    X_validation_pfi = X_validation_all[[col for col in pfi_features if col in X_validation_all.columns]]
    X_validation_lasso = X_validation_all[[col for col in lasso_features if col in X_validation_all.columns]]
    sample_ids_bioproject = group[['Run_ID', 'BioProject', 'Phenotype']]

    y_validation_rfe = y_validation.loc[X_validation_rfe.index]
    y_validation_pfi = y_validation.loc[X_validation_pfi.index]
    y_validation_lasso = y_validation.loc[X_validation_lasso.index]

    print(f"BioProject {project_id} - y_validation unique values: {np.unique(y_validation)}")
    print(f"BioProject {project_id} - X_validation shape: {X_validation.shape}")
    print(f"BioProject {project_id} - X_validation_rfe shape: {X_validation_rfe.shape}")
    print(f"BioProject {project_id} - y_validation_rfe shape: {y_validation_rfe.shape}")

    for models, X_set, y_set, method in [
        (models_all_features, X_validation, y_validation, "All_Features"),
        (models_pfi, X_validation_pfi, y_validation_pfi, "PFI"),
        (models_rfe, X_validation_rfe, y_validation_rfe, "RFE"),
        (models_lasso, X_validation_lasso, y_validation_lasso, "LASSO")
    ]:
        for model_name, model in models.items():
            matrix = evaluate_model(model, X_set, y_set, model.predict_proba(X_set)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_set))
            matrix.update({'BioProject': project_id, 'Phenotype': phenotype_str, 'Model': model_name,
                            'Feature_Selection': method, 'BioProject_Phenotype_Counts': phenotype_counts_str})
            all_metrics.append(matrix)

        all_bioproject_feature_importance.extend(evaluate_feature_importance(models, X_set, y_set, method, "bioproject", featrue_importance_dir, project_id))
        all_cm_results.extend(compute_confusion_matrix(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence, project_id))
        all_classification_results.extend(sample_classification_details(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence, project_id))
        all_cm_results.extend(plot_confusion_matrix(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence, project_id))
        all_misclass_results.extend(analyze_misclassifications(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence))
        all_classification_reports.extend(generate_classification_report(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence, project_id))
        all_classification_reports.extend(plot_classification_report(models, X_set, y_set, sample_ids_bioproject, method, "bioproject", output_dir, prevalence, project_id))

metrics_df = pd.DataFrame(all_metrics)
metrics_df.to_csv(os.path.join(output_dir, f'model_evaluation_on_each_bioproject_validation_datasets_{prevalence}.csv'), index=False)
print(f"BioProject evaluation results saved to: {os.path.join(output_dir, f'model_evaluation_on_each_bioproject_validation_datasets_{prevalence}.csv')}")

if all_bioproject_feature_importance:
    pd.concat(all_bioproject_feature_importance, ignore_index=True).to_csv(os.path.join(output_dir, f'feature_importance_bioproject_{prevalence}.tsv'), sep='\t', index=False)
pd.concat(all_cm_results, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'confusion_matrix_bioproject_{prevalence}.tsv'), sep='\t', index=False)
pd.concat(all_classification_results, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'classification_details_bioproject_{prevalence}.tsv'), sep='\t', index=False)
if all_misclass_results:
    pd.concat(all_misclass_results, ignore_index=True).to_csv(os.path.join(output_dir, 'confusion_matrix_results', f'misclassifications_bioproject_{prevalence}.tsv'), sep='\t', index=False)
if all_classification_reports:
    pd.concat(all_classification_reports, ignore_index=True).to_csv(os.path.join(output_dir, 'classification_reports', f'classification_reports_bioproject_{prevalence}.tsv'), sep='\t', index=False)

for project_id, group in grouped_data:
    phenotype_counts = group["Phenotype"].value_counts().to_dict()
    phenotype_counts_str = f"{project_id}: {','.join([f'{key}({value})' for key, value in phenotype_counts.items()])}"
    bioproject_to_phenotype[project_id] = phenotype_counts_str

    validation_df = group.drop(columns=[col for col in columns_to_drop if col in group.columns])
    y_validation = validation_df['Class_Label'].astype(int)
    X_validation_all = validation_df.drop(columns=['Class_Label'])

    X_validation = X_validation_all[[col for col in X_train.columns if col in X_validation_all.columns]]
    X_validation_rfe = X_validation_all[[col for col in rfe_features if col in X_validation_all.columns]].dropna()
    X_validation_pfi = X_validation_all[[col for col in pfi_features if col in X_validation_all.columns]]
    X_validation_lasso = X_validation_all[[col for col in lasso_features if col in X_validation_all.columns]]

    y_validation_rfe = y_validation.loc[X_validation_rfe.index]
    y_validation_pfi = y_validation.loc[X_validation_pfi.index]
    y_validation_lasso = y_validation.loc[X_validation_lasso.index]

    for models, X_set, y_set in [
        (models_all_features, X_validation, y_validation),
        (models_rfe, X_validation_rfe, y_validation_rfe),
        (models_pfi, X_validation_pfi, y_validation_pfi),
        (models_lasso, X_validation_lasso, y_validation_lasso)
    ]:
        roc_data = plot_roc_curves_bioproject(models, X_set, y_set, project_id, phenotype_counts_str, roc_data_dir)
        if not roc_data.empty:
            all_roc_data.append(roc_data)

if all_roc_data:
    all_roc_data_df = pd.concat(all_roc_data, ignore_index=True)
    all_roc_data_df.to_csv(os.path.join(roc_data_dir, f'roc_data_bioproject_{prevalence}.tsv'), sep='\t', index=False)

    g = sns.FacetGrid(all_roc_data_df, col="BioProject", hue="Model", col_wrap=4, height=6, aspect=1.5)
    g.map(sns.lineplot, "FPR", "TPR", marker="o")
    g.add_legend()
    g.set_axis_labels("False Positive Rate", "True Positive Rate")
    for ax, project_id in zip(g.axes.flat, g.col_names):
        ax.set_title(bioproject_to_phenotype.get(project_id, project_id), fontsize=12, pad=10)
    g.fig.tight_layout(pad=1.5)
    output_path = os.path.join(roc_data_dir, f'roc_facet_grid_validation_{prevalence}.png')
    g.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"ROC facet grid saved to {output_path}")
    plt.close(g.fig)
else:
    print("No ROC data generated due to single-class BioProjects")

print(f"\nEvaluation pipeline complete for prevalence: {prevalence}")
print(f"All outputs written under: {output_dir}")