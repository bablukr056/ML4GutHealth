#!/lustrehome/babluuniba2022/miniconda3/envs/ml/bin/python

# ==============================================================================
# Logistic Regression (ElasticNet) classifier training on prevalence-filtered,
# CLR-transformed species-level abundance data, comparing four
# feature-selection strategies: all-features, RFECV, Lasso, and permutation
# importance.
#
# The input filename (a prevalence-threshold-specific merged OTU/metadata
# file) is supplied as a command-line argument at runtime. Model training,
# including hyperparameter search (RandomizedSearchCV followed by
# GridSearchCV), was run as an HTCondor job on the ReCaS-Bari HPC cluster,
# with 32 CPU cores and 272 GB of RAM allocated per job.
# ==============================================================================

import argparse
import sklearn
import numpy as np
import pandas as pd
import joblib
import os
import time as time_module
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import accuracy_score, classification_report, make_scorer
from dask_ml.model_selection import RandomizedSearchCV, GridSearchCV
from sklearn.metrics import f1_score

print("NumPy version:", np.__version__)
print("scikit-learn version:", sklearn.__version__)
print("Pandas version:", pd.__version__)


def logistic_regression_model(X_train, X_test, y_train, y_test, model_name, output_dir, model_dir):
    try:
        model = LogisticRegression(solver='saga', penalty='elasticnet', max_iter=10000, class_weight="balanced")
        KFCV = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        f1 = make_scorer(f1_score, average='binary')

        param_dist = {
            'C': np.logspace(-4, 4, 20),
            'l1_ratio': np.linspace(0, 1, 11)  # ElasticNet: l1_ratio between 0 (Ridge) and 1 (Lasso)
        }

        random_search = RandomizedSearchCV(model, param_distributions=param_dist, n_iter=60, cv=KFCV,
                                           scoring=f1, n_jobs=-1, random_state=42)
        start = time_module.time()
        random_search.fit(X_train, y_train)
        print(f"Random search for {model_name} took {time_module.time() - start:.2f} seconds")
        print(f"Best parameters from random search: {random_search.best_params_}")
        print(f"Best F1 score from random search: {random_search.best_score_:.4f}")

        best_params = random_search.best_params_
        pd.DataFrame(random_search.cv_results_).to_csv(os.path.join(model_dir, f"LOGREG_RANDOM_SEARCH_{model_name}.tsv"), sep="\t", index=False)

        param_grid = {
            'C': [best_params['C'] * 0.1, best_params['C'], best_params['C'] * 2],
            'l1_ratio': [best_params['l1_ratio']]
        }

        grid_search = GridSearchCV(model, param_grid, cv=KFCV, scoring=f1, n_jobs=-1)
        start = time_module.time()
        grid_search.fit(X_train, y_train)
        print(f"Grid search for {model_name} took {time_module.time() - start:.2f} seconds")
        pd.DataFrame(grid_search.cv_results_).to_csv(os.path.join(model_dir, f"LOGREG_GRID_SEARCH_{model_name}.tsv"), sep="\t", index=False)
        best_model = grid_search.best_estimator_

        output_model_path = os.path.join(model_dir, f"LOGREG_{model_name}.joblib")
        joblib.dump(best_model, output_model_path)

        predictions = best_model.predict(X_test)
        accuracy = accuracy_score(y_test, predictions)
        print(f"Test accuracy for {model_name}: {accuracy:.4f}")
        print(f"Classification report for {model_name}:")
        print(classification_report(y_test, predictions))

        return best_model
    except Exception as e:
        print(f"Error in logistic_regression_model for {model_name}: {str(e)}")
        return None


def process_prevalence_file(filename):
    try:
        base_dir = "/lustrehome/babluuniba2022/bablu/condor/Data_Analysis/JUNE_2025/prevalence_based_threshold"
        base_model_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/"

        logreg_model_dir = os.path.join(base_model_dir, "logreg_models")
        rfecv_feature_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/svm_models_full"
        lasso_feature_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/lasso_models"
        pfi_feature_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/"

        os.makedirs(logreg_model_dir, exist_ok=True)

        data_path = os.path.join(base_dir, filename)
        prevalence_file = os.path.splitext(os.path.basename(data_path))[0]
        prevalence_suffix = prevalence_file.replace("clr_otu_metadata_merged_prevalence_bacterial_", "").replace("clr_otu_metadata_merged_prevalence_", "")
        print(f"\nProcessing prevalence file: {prevalence_file} with suffix: {prevalence_suffix}")

        df = pd.read_csv(data_path, sep='\t', low_memory=False)

        columns_to_drop = ["Run_ID", "BioProject", "BioSample", "Health_status", "Phenotype",
                           "Full_Name", "Sex", "Age", "Location", "Sample", "BMI", "Platform", "Author"]
        columns_to_drop = [col for col in columns_to_drop if col in df.columns]
        df = df.drop(columns=columns_to_drop)

        if 'Class_Label' not in df.columns:
            raise ValueError("Class_Label column not found in the dataset")

        y = df['Class_Label']
        X = df.drop(columns=['Class_Label'])

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

        model_name_prefix = f"BACTERIAL_{prevalence_suffix}"

        # --------------------------------------------------------------------
        # Logistic Regression (ElasticNet) trained on the full (all-species)
        # feature set.
        # --------------------------------------------------------------------
        print("\nTraining Logistic Regression on all features...")
        logistic_regression_model(X_train, X_test, y_train, y_test,
                                  model_name=f"{model_name_prefix}_ALL_FEATURES",
                                  output_dir=base_model_dir,
                                  model_dir=logreg_model_dir)

        # --------------------------------------------------------------------
        # RFECV-selected features (pre-computed and loaded from a separate
        # SVM-based RFECV run; see svm_model_training.py). RFECV itself was
        # fit on the training set only, to avoid data leakage.
        # --------------------------------------------------------------------
        rfecv_feature_file = os.path.join(rfecv_feature_dir, f"RFECV_SELECTED_FEATURES_WITH_COEFFICIENTS_{prevalence_suffix}.tsv")
        if os.path.exists(rfecv_feature_file):
            rfecv_df = pd.read_csv(rfecv_feature_file, sep='\t')
            selected_features = rfecv_df['Feature_Name'].tolist()
            selected_features = [f for f in selected_features if f in X.columns]
            X_train_rfecv = X_train[selected_features]
            X_test_rfecv = X_test[selected_features]

            print(f"\nTraining Logistic Regression on RFECV-selected features ({len(selected_features)} features)...")
            logistic_regression_model(X_train_rfecv, X_test_rfecv, y_train, y_test,
                                      model_name=f"{model_name_prefix}_RFECV_FEATURES",
                                      output_dir=base_model_dir,
                                      model_dir=logreg_model_dir)
        else:
            print("\nRFECV feature file not found, skipping RFECV-based Logistic Regression training.")

        # --------------------------------------------------------------------
        # Lasso-selected features (pre-computed and loaded from a separate
        # Lasso run).
        # --------------------------------------------------------------------
        lasso_feature_file = os.path.join(lasso_feature_dir, f"lasso_based_featrues_on_{prevalence_suffix}.tsv")
        if os.path.exists(lasso_feature_file):
            lasso_df = pd.read_csv(lasso_feature_file, sep='\t')
            selected_features = lasso_df['Feature'].tolist()
            selected_features = [f for f in selected_features if f in X.columns]
            X_train_lasso = X_train[selected_features]
            X_test_lasso = X_test[selected_features]

            print(f"\nTraining Logistic Regression on LASSO-selected features ({len(selected_features)} features)...")
            logistic_regression_model(X_train_lasso, X_test_lasso, y_train, y_test,
                                      model_name=f"{model_name_prefix}_LASSO_FEATURES",
                                      output_dir=base_model_dir,
                                      model_dir=logreg_model_dir)
        else:
            print("\nLASSO feature file not found, skipping LASSO-based Logistic Regression training.")

        # --------------------------------------------------------------------
        # Permutation Importance-selected features. Importance scores were
        # computed on the held-out test set of a previously trained Random
        # Forest model (see rf_model_training.py).
        # --------------------------------------------------------------------
        pfi_feature_file = os.path.join(pfi_feature_dir, f"PERMUTATION_IMPORTANCE_BACTERIAL_bacterial_{prevalence_suffix}.tsv")
        if os.path.exists(pfi_feature_file):
            pfi_df = pd.read_csv(pfi_feature_file, sep='\t')
            selected_features = pfi_df[pfi_df['Importance_Mean'] > 0]['Feature'].tolist()
            selected_features = [f for f in selected_features if f in X.columns]
            X_train_pfi = X_train[selected_features]
            X_test_pfi = X_test[selected_features]

            print(f"\nTraining Logistic Regression on Permutation Importance-selected features ({len(selected_features)} features)...")
            logistic_regression_model(X_train_pfi, X_test_pfi, y_train, y_test,
                                      model_name=f"{model_name_prefix}_PERMUTATION_FEATURES",
                                      output_dir=base_model_dir,
                                      model_dir=logreg_model_dir)
        else:
            print("\nPermutation Importance feature file not found, skipping PFI-based Logistic Regression training.")

    except Exception as e:
        print(f"Error processing {filename}: {str(e)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Train Logistic Regression model (ElasticNet) on multiple feature "
                     "selection strategies (all-features, RFECV, Lasso, Permutation Importance). "
                     "The input filename is supplied via the terminal at runtime; model ")
    parser.add_argument("filename", help="Input data file inside prevalence_based_threshold directory,  from the terminal at runtime.")
    args = parser.parse_args()

    process_prevalence_file(args.filename)
