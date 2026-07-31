#! /lustrehome/babluuniba2022/miniconda3/envs/ml/bin/python

# ==============================================================================
# SVM classifier training (RBF and Linear kernels) on prevalence-filtered,
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
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import accuracy_score, classification_report, make_scorer, matthews_corrcoef
from sklearn.feature_selection import RFECV
from dask_ml.model_selection import RandomizedSearchCV, GridSearchCV

print("NumPy version:", np.__version__)
print("scikit-learn version:", sklearn.__version__)
print("Pandas version:", pd.__version__)


def svm_model_building(X_train, X_test, y_train, y_test, model_name, kernel, output_dir, model_dir):
    try:
        svm = SVC(class_weight="balanced")
        KFCV = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        mcc = make_scorer(matthews_corrcoef)

        if kernel == 'rbf':
            param_dist_random = {
                'C': [0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 4, 10],
                'gamma': [1e-4, 1e-3, 1e-2, 0.1, 0.5, 1, 'scale', 'auto'],
                'kernel': ['rbf']
            }
        elif kernel == 'linear':
            param_dist_random = {
                'C': [0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 4, 10],
                'kernel': ['linear']
            }
        else:
            raise ValueError(f"Unsupported kernel: {kernel}")

        random_search = RandomizedSearchCV(svm, param_distributions=param_dist_random, n_iter=50, cv=KFCV, scoring=mcc, n_jobs=-1, random_state=42)
        start = time_module.time()
        random_search.fit(X_train, y_train)
        print(f"Random search for {model_name} ({kernel} kernel) took {time_module.time() - start:.2f} seconds")

        random_search_results = pd.DataFrame(random_search.cv_results_)
        random_search_output = os.path.join(output_dir, f"parameters_randomsearchcv_{model_name}_{kernel}.tsv")
        random_search_results.to_csv(random_search_output, sep="\t", index=False)

        best_params_random = random_search.best_params_

        if kernel == 'rbf':
            if isinstance(best_params_random['gamma'], (int, float)):
                gamma_grid = [best_params_random['gamma'] * 0.1, best_params_random['gamma'], best_params_random['gamma'] * 10]
            else:
                gamma_grid = [best_params_random['gamma']]
                print(f"Best gamma is {best_params_random['gamma']}, using it directly in GridSearchCV")
            param_grid = {
                'C': [best_params_random['C'] * 0.1, best_params_random['C'], best_params_random['C'] * 2],
                'gamma': gamma_grid,
                'kernel': ['rbf']
            }
        else:
            param_grid = {
                'C': [best_params_random['C'] * 0.1, best_params_random['C'], best_params_random['C'] * 2],
                'kernel': ['linear']
            }

        grid_search = GridSearchCV(svm, param_grid, cv=KFCV, scoring=mcc, n_jobs=-1)
        start = time_module.time()
        grid_search.fit(X_train, y_train)
        print(f"Grid search for {model_name} ({kernel} kernel) took {time_module.time() - start:.2f} seconds")

        grid_search_results = pd.DataFrame(grid_search.cv_results_)
        grid_search_output = os.path.join(output_dir, f"parameters_gridsearchcv_{model_name}_{kernel}.tsv")
        grid_search_results.to_csv(grid_search_output, sep="\t", index=False)

        best_model = grid_search.best_estimator_
        model_output_file = os.path.join(model_dir, f"SVM_{model_name}_{kernel}.joblib")
        joblib.dump(best_model, model_output_file)

        predictions = best_model.predict(X_test)
        test_accuracy = accuracy_score(y_test, predictions)
        print(f"Test accuracy for {model_name} ({kernel} kernel): {test_accuracy:.4f}")
        print(f"Classification report for {model_name} ({kernel} kernel):\n{classification_report(y_test, predictions)}")

        return best_model
    except Exception as e:
        print(f"Error in svm_model_building for {model_name} ({kernel} kernel): {str(e)}")
        return None


def process_prevalence_file(filename):
    try:
        base_dir = "/lustrehome/babluuniba2022/bablu/condor/Data_Analysis/JUNE_2025/prevalence_based_threshold"
        data_path = os.path.join(base_dir, filename)
        output_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev"
        model_dir = os.path.join(output_dir, "svm_models_full")
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(model_dir, exist_ok=True)

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

        # Single train-test split for all methods
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

        model_name_prefix = f"BACTERIAL_{prevalence_suffix}"

        # --------------------------------------------------------------------
        # SVM trained on the full (all-species) feature set.
        # --------------------------------------------------------------------
        print("\nTraining SVM on all features...")
        for kernel in ['rbf', 'linear']:
            svm_model_building(X_train, X_test, y_train, y_test,
                               model_name=f"{model_name_prefix}_ALL_FEATURES",
                               kernel=kernel,
                               output_dir=output_dir,
                               model_dir=model_dir)

        # --------------------------------------------------------------------
        # RFECV feature selection, based on the linear SVM model above.
        # Fit on the training set only, to avoid data leakage into the
        # held-out test set.
        # --------------------------------------------------------------------
        linear_model_path = os.path.join(model_dir, f"SVM_{model_name_prefix}_ALL_FEATURES_linear.joblib")
        if not os.path.exists(linear_model_path):
            print(f"Linear SVM model file not found for RFECV: {linear_model_path}")
        else:
            print("\nStarting RFECV feature selection using linear SVM model...")
            SVM_LINEAR = joblib.load(linear_model_path)

            mcc_scorer = make_scorer(matthews_corrcoef)
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

            start_time = time_module.time()
            rfecv = RFECV(estimator=SVM_LINEAR, step=1, cv=cv, scoring=mcc_scorer, n_jobs=-1, verbose=1)
            rfecv.fit(X_train, y_train)  # Use train data for feature selection
            end_time = time_module.time()
            print(f"Time taken for RFECV: {end_time - start_time:.2f} seconds")

            print("Optimal number of features:", rfecv.n_features_)

            best_features = X_train.columns[rfecv.support_]
            coefficients = SVM_LINEAR.coef_[0][rfecv.support_]
            ref_cof = pd.DataFrame({'Feature_Name': best_features, 'Coefficient': coefficients})
            ref_cof_path = os.path.join(output_dir, f"RFECV_SELECTED_FEATURES_WITH_COEFFICIENTS_{prevalence_suffix}.tsv")
            ref_cof.to_csv(ref_cof_path, sep='\t', index=False)
            print(f"RFECV selected features with coefficients saved to: {ref_cof_path}")

            X_train_rfecv = X_train[best_features]
            X_test_rfecv = X_test[best_features]

            print("\nTraining SVM on RFECV selected features...")
            for kernel in ['rbf', 'linear']:
                svm_model_building(X_train_rfecv, X_test_rfecv, y_train, y_test,
                                   model_name=f"{model_name_prefix}_RFECV_FEATURES",
                                   kernel=kernel,
                                   output_dir=output_dir,
                                   model_dir=model_dir)

        # --------------------------------------------------------------------
        # Lasso-based SVM training
        # --------------------------------------------------------------------
        lasso_features_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/lasso_models"
        lasso_feature_file = f"lasso_based_featrues_on_{prevalence_suffix}.tsv"
        lasso_feature_path = os.path.join(lasso_features_dir, lasso_feature_file)

        if os.path.exists(lasso_feature_path):
            print(f"\nFound Lasso selected features file: {lasso_feature_path}")
            lasso_features_df = pd.read_csv(lasso_feature_path, sep='\t')
            selected_features = lasso_features_df['Feature'].tolist()
            selected_features = [f for f in selected_features if f in X.columns]
            if not selected_features:
                print("No valid Lasso-selected features found in dataset columns, skipping Lasso-based SVM training.")
            else:
                X_train_lasso = X_train[selected_features]
                X_test_lasso = X_test[selected_features]

                print(f"\nTraining SVM on Lasso-selected features ({len(selected_features)} features)...")
                for kernel in ['rbf', 'linear']:
                    svm_model_building(X_train_lasso, X_test_lasso, y_train, y_test,
                                       model_name=f"{model_name_prefix}_LASSO_FEATURES",
                                       kernel=kernel,
                                       output_dir=output_dir,
                                       model_dir=model_dir)
        else:
            print(f"\nLasso selected features file not found for prevalence suffix '{prevalence_suffix}'. Skipping Lasso-based SVM training.")

        # --------------------------------------------------------------------
        # Permutation Importance-based SVM training. Uses feature-importance
        # scores computed on the held-out test set of a previously trained
        # Random Forest model (see rf_model_training.py).
        # --------------------------------------------------------------------
        perm_features_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev"
        perm_feature_file = f"PERMUTATION_IMPORTANCE_BACTERIAL_bacterial_{prevalence_suffix}.tsv"
        perm_feature_path = os.path.join(perm_features_dir, perm_feature_file)

        if os.path.exists(perm_feature_path):
            print(f"\nFound Permutation Importance selected features file: {perm_feature_path}")
            perm_features_df = pd.read_csv(perm_feature_path, sep='\t')
            selected_pfi = perm_features_df[perm_features_df['Importance_Mean'] > 0]['Feature'].tolist()
            selected_pfi = [f for f in selected_pfi if f in X.columns]
            if not selected_pfi:
                print("No valid Permutation Importance-selected features with Importance_Mean > 0 found in dataset columns, skipping Permutation-based SVM training.")
            else:
                X_train_perm = X_train[selected_pfi]
                X_test_perm = X_test[selected_pfi]

                print(f"\nTraining SVM on Permutation Importance-selected features ({len(selected_pfi)} features)...")
                for kernel in ['rbf', 'linear']:
                    svm_model_building(X_train_perm, X_test_perm, y_train, y_test,
                                       model_name=f"{model_name_prefix}_PERMUTATION_FEATURES",
                                       kernel=kernel,
                                       output_dir=output_dir,
                                       model_dir=model_dir)
        else:
            print(f"\nPermutation Importance selected features file not found for prevalence suffix '{prevalence_suffix}'. Skipping Permutation-based SVM training.")

    except Exception as e:
        print(f"Error processing {filename}: {str(e)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run SVM classification pipeline with RBF and Linear kernels  ")
    parser.add_argument("filename", help="Filename inside the prevalence_based_threshold directory, supplied from the terminal at runtime.")
    args = parser.parse_args()

    process_prevalence_file(args.filename)
