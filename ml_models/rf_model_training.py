#! /lustrehome/babluuniba2022/miniconda3/envs/ml/bin/python

# ==============================================================================
# Random Forest classifier training on prevalence-filtered, CLR-transformed
# species-level abundance data, comparing four feature-selection strategies:
# all-features, permutation importance, RFECV, and Lasso.
#
# The input filename (a prevalence-threshold-specific merged OTU/metadata
# python rf_model_training.py clr_otu_metadata_merged_prevalence_bacterial_20pct.tsv
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
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV, RandomizedSearchCV
from sklearn.metrics import accuracy_score, classification_report, make_scorer, matthews_corrcoef
from sklearn.inspection import permutation_importance

print("NumPy version:", np.__version__)
print("scikit-learn version:", sklearn.__version__)
print("Pandas version:", pd.__version__)


def rf_model_building(X_train, X_test, y_train, y_test, model_name, output_dir, model_dir):
    try:
        rf = RandomForestClassifier(random_state=42, n_jobs=-1, oob_score=True, class_weight="balanced")
        KFCV = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        mcc = make_scorer(matthews_corrcoef)

        random_search_grid = {
            "n_estimators": list(range(1000, 1200, 50)),
            "max_features": list(range(10, 100, 20)),
            "min_samples_leaf": list(range(10, 100, 10)),
            "max_depth": [None, 5, 10, 15],
            "class_weight": ['balanced']
        }

        random_search = RandomizedSearchCV(estimator=rf, param_distributions=random_search_grid, n_iter=50,
                                          cv=KFCV, scoring=mcc, n_jobs=-1, random_state=42, return_train_score=True)
        start = time_module.time()
        random_search.fit(X_train, y_train)
        print(f"Random search for {model_name} took {time_module.time() - start:.2f} seconds")
        print(f"Best parameters from random search: {random_search.best_params_}")
        print(f"Best MCC score from random search: {random_search.best_score_:.4f}")

        random_search_results = pd.DataFrame(random_search.cv_results_)
        random_search_output = os.path.join(output_dir, f"parameters_randomsearchcv_{model_name}.tsv")
        random_search_results.to_csv(random_search_output, sep="\t", index=False)

        best_params_random = random_search.best_params_
        param_grid = {
            "n_estimators": list(range(max(100, best_params_random["n_estimators"] - 50), best_params_random["n_estimators"] + 50, 30)),
            "max_features": list(range(max(1, best_params_random["max_features"] - 5), best_params_random["max_features"] + 5, 3)),
            "min_samples_leaf": list(range(max(1, best_params_random["min_samples_leaf"] - 5), best_params_random["min_samples_leaf"] + 5, 4)),
            "max_depth": [best_params_random["max_depth"]] if best_params_random["max_depth"] else [None],
            "class_weight": [best_params_random["class_weight"]]
        }

        grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=KFCV, scoring=mcc, n_jobs=-1)
        start = time_module.time()
        grid_search.fit(X_train, y_train)
        print(f"Grid search for {model_name} took {time_module.time() - start:.2f} seconds")

        grid_search_results = pd.DataFrame(grid_search.cv_results_)
        grid_search_output = os.path.join(output_dir, f"parameters_gridsearchcv_{model_name}.tsv")
        grid_search_results.to_csv(grid_search_output, sep="\t", index=False)

        best_model = grid_search.best_estimator_
        model_output_file = os.path.join(model_dir, f"RF_{model_name}.joblib")
        joblib.dump(best_model, model_output_file)

        predictions = best_model.predict(X_test)
        test_accuracy = accuracy_score(y_test, predictions)
        print(f"Test accuracy for {model_name}: {test_accuracy:.4f}")
        print(f"Classification report for {model_name}:\n{classification_report(y_test, predictions)}")

        return best_model
    except Exception as e:
        print(f"Error in rf_model_building for {model_name}: {str(e)}")
        return None


def process_prevalence_file(filename):
    try:
        base_dir = "/lustrehome/babluuniba2022/bablu/condor/Data_Analysis/JUNE_2025/prevalence_based_threshold"
        output_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev"
        model_dir = os.path.join(output_dir, "random_forest_models")
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(model_dir, exist_ok=True)

        prevalence_file = os.path.splitext(os.path.basename(filename))[0]
        prevalence_suffix = prevalence_file.replace("clr_otu_metadata_merged_prevalence_bacterial_", "")
        print(f"\nProcessing prevalence file: {prevalence_file} with suffix: {prevalence_suffix}")

        data_path = os.path.join(base_dir, filename)
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
        print(f"Training set size: {X_train.shape[0]} samples, Test set size: {X_test.shape[0]} samples")

        model_name_prefix = f"BACTERIAL_{prevalence_suffix}"

        # Train RF on all features
        print("\nTraining RF on all features...")
        base_model = rf_model_building(X_train, X_test, y_train, y_test,
                                       model_name=f"{model_name_prefix}_ALL_FEATURES",
                                       output_dir=output_dir, model_dir=model_dir)

        # --- PERMUTATION IMPORTANCE-BASED RF TRAINING ---   # Ref: https://scikit-learn.org/stable/modules/permutation_importance.html
        if base_model is not None:
            model_path = os.path.join(model_dir, f"RF_{model_name_prefix}_ALL_FEATURES.joblib")
            if not os.path.exists(model_path):
                print(f"RF model file not found for PFI: {model_path}")
            else:
                print(f"\nPerforming permutation importance for {model_name_prefix}...")
                rf_model = joblib.load(model_path)
                result = permutation_importance(rf_model, X_test, y_test, n_repeats=10, random_state=42, n_jobs=-1)
                importance_df = pd.DataFrame({
                    'Feature': X_train.columns,
                    'Importance_Mean': result.importances_mean,
                    'Importance_Std': result.importances_std
                })

                selected_features = importance_df[importance_df['Importance_Mean'] > 0]['Feature'].tolist()
                perm_output = os.path.join(output_dir, f"PERMUTATION_IMPORTANCE_{model_name_prefix}.tsv")
                importance_df.to_csv(perm_output, sep='\t', index=False)

                print(f"Selected {len(selected_features)} features with importance > 0")
                if selected_features:
                    X_train_perm = X_train[selected_features]
                    X_test_perm = X_test[selected_features]

                    print(f"\nTraining RF on Permutation Importance-selected features ({len(selected_features)} features)...")
                    rf_model_building(X_train_perm, X_test_perm, y_train, y_test,
                                      model_name=f"{model_name_prefix}_PERMUTATION_FEATURES",
                                      output_dir=output_dir, model_dir=model_dir)
                else:
                    print("No valid Permutation Importance-selected features with Importance_Mean > 0 found, skipping PFI-based RF training.")

        # --- RFECV FEATURE SELECTION (LOAD PRE-SELECTED FEATURES) ---
        rfecv_features_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/svm_models_full"
        rfecv_feature_file = f"RFECV_SELECTED_FEATURES_WITH_COEFFICIENTS_{prevalence_suffix}.tsv"
        rfecv_feature_path = os.path.join(rfecv_features_dir, rfecv_feature_file)

        if os.path.exists(rfecv_feature_path):
            print(f"\nFound RFECV selected features file: {rfecv_feature_path}")
            rfecv_features_df = pd.read_csv(rfecv_feature_path, sep='\t')
            selected_features = rfecv_features_df['Feature_Name'].tolist()
            selected_features = [f for f in selected_features if f in X.columns]
            if not selected_features:
                print("No valid RFECV-selected features found in dataset columns, skipping RFECV-based RF training.")
            else:
                X_train_rfecv = X_train[selected_features]
                X_test_rfecv = X_test[selected_features]

                print(f"\nTraining RF on RFECV-selected features ({len(selected_features)} features)...")
                rf_model_building(X_train_rfecv, X_test_rfecv, y_train, y_test,
                                  model_name=f"{model_name_prefix}_RFECV_FEATURES",
                                  output_dir=output_dir, model_dir=model_dir)
        else:
            print(f"\nRFECV selected features file not found for prevalence suffix '{prevalence_suffix}'. Skipping RFECV-based RF training.")

        # --- LASSO-BASED RF TRAINING ---
        lasso_features_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev/lasso_models"
        lasso_feature_file = f"lasso_based_featrues_on_{prevalence_suffix}.tsv"
        lasso_feature_path = os.path.join(lasso_features_dir, lasso_feature_file)

        if os.path.exists(lasso_feature_path):
            print(f"\nFound Lasso selected features file: {lasso_feature_path}")
            lasso_features_df = pd.read_csv(lasso_feature_path, sep='\t')
            selected_features = lasso_features_df['Feature'].tolist()
            selected_features = [f for f in selected_features if f in X.columns]
            if not selected_features:
                print("No valid Lasso-selected features found in dataset columns, skipping Lasso-based RF training.")
            else:
                X_train_lasso = X_train[selected_features]
                X_test_lasso = X_test[selected_features]

                print(f"\nTraining RF on Lasso-selected features ({len(selected_features)} features)...")
                rf_model_building(X_train_lasso, X_test_lasso, y_train, y_test,
                                  model_name=f"{model_name_prefix}_LASSO_FEATURES",
                                  output_dir=output_dir, model_dir=model_dir)
        else:
            print(f"\nLasso selected features file not found for prevalence suffix '{prevalence_suffix}'. Skipping Lasso-based RF training.")

    except Exception as e:
        print(f"Error processing {filename}: {str(e)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run RF classification pipeline on prevalence-filtered bacterial OTU data, "
                     "including all-features, RFECV, Lasso, and Permutation Importance feature "
                     "selection. The input filename is supplied via the terminal at runtime; "
                     "model training was performed as an HTCondor job with 32 CPU cores and "
                     "272 GB RAM."
    )
    parser.add_argument("filename", help="Filename inside the prevalence_based_threshold directory, supplied from the terminal at runtime.")
    args = parser.parse_args()

    process_prevalence_file(args.filename)
