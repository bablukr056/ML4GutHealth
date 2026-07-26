#! /lustrehome/babluuniba2022/miniconda3/envs/ml/bin/python


# ==============================================================================
#
# This notebook evaluates whether normalization strategy (Raw counts, TSS, CLR) and CLR pseudocount value (1.0, 0.5, 0.1) influence the predictive performance of our best-performing model (SVM-RBF, PFI-selected features), using the full hyperparameter-tuning pipeline: RandomizedSearchCV (50 iterations, 5-fold stratified CV, F1 scoring) followed by a refined GridSearchCV around the best parameters found.
#
# Five configurations are evaluated:
# 1. Raw counts (no normalization)
# 2. TSS (Total Sum Scaling / relative abundance)
# 3. CLR, pseudocount = 1.0 (value used in the primary manuscript analysis)
# 4. CLR, pseudocount = 0.5
# 5. CLR, pseudocount = 0.1
#
# All input data is loaded once at the start. All outputs (per-configuration tuning-result files, summary table, and one figure per technique) are saved to a single, clearly named output folder.
# ==============================================================================

import os
import time as time_module
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.model_selection import train_test_split, StratifiedKFold, RandomizedSearchCV, GridSearchCV
from sklearn.svm import SVC
from sklearn.metrics import f1_score, roc_auc_score

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.size'] = 11

OUTPUT_DIR = "result_with_normalization" # change the output directory name as needed
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==============================================================================
#  Section 1 — Load Input Data (loaded once, reused throughout)
# ==============================================================================

meta_data = pd.read_csv("otu_data/metadata_for_7452.tsv", sep="\t")

training_corrected = pd.read_csv("otu_data/corrected_training_7452.csv", sep="\t", header=0, index_col=0, engine="c")

pfi_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev"
pfi_file = os.path.join(pfi_dir, "PERMUTATION_IMPORTANCE_BACTERIAL_bacterial_20pct.tsv")
pfi_results = pd.read_csv(pfi_file, sep="\t")

SPECIES_COLUMN = pfi_results.columns[0]
IMPORTANCE_COLUMN = pfi_results.columns[1]

pfi_species_list = pfi_results[pfi_results[IMPORTANCE_COLUMN] > 0][SPECIES_COLUMN].tolist()

# ==============================================================================
#  Section 2 — Prepare Class Labels
# ==============================================================================

meta_data_indexed = meta_data.set_index("Run_ID")
sample_ids = training_corrected.columns.tolist()
health_labels = meta_data_indexed.loc[sample_ids, "Health_status"]

label_map = {"HEALTHY": 1, "NON-HEALTHY": 0} # label mapping for out target variable
y = health_labels.map(label_map).values

# ==============================================================================
#  Section 3 — Define Normalization / CLR-Pseudocount Configurations
# ==============================================================================

def normalize_raw_counts(count_df):
    """CONFIGURATION 1: Raw counts (no normalization applied)."""
    return count_df.copy().astype(float)


def normalize_tss(count_df):
    """CONFIGURATION 2: Total Sum Scaling (relative abundance) - each
    sample's counts divided by that sample's total read count."""
    df = count_df.copy().astype(float)
    sample_totals = df.sum(axis=0)
    return df.divide(sample_totals, axis=1)


def normalize_clr(count_df, pseudocount):
    """CONFIGURATIONS 3-5: Centered Log-Ratio transformation with an
    additive pseudocount for zero-value handling (manuscript method),
    evaluated at three pseudocount values."""
    df = count_df.copy().astype(float)
    df = df + pseudocount
    log_df = np.log(df)
    geometric_mean = log_df.mean(axis=0)
    return log_df.subtract(geometric_mean, axis=1)


configurations = {
    "Raw_Counts": normalize_raw_counts(training_corrected),
    "TSS_Relative_Abundance": normalize_tss(training_corrected),
    "CLR_pseudocount_1.0": normalize_clr(training_corrected, 1.0),
    "CLR_pseudocount_0.5": normalize_clr(training_corrected, 0.5),
    "CLR_pseudocount_0.1": normalize_clr(training_corrected, 0.1),
}

# ==============================================================================
#  Section 4 — Hyperparameter-Tuning Function
#
# ==============================================================================

def tune_and_evaluate_svm_rbf(X_train, X_test, y_train, y_test, config_name, output_dir):
    svm = SVC(class_weight="balanced")
    KFCV = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    param_dist_random = {
        'C': [0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 4, 10],
        'gamma': [1e-4, 1e-3, 1e-2, 0.1, 0.5, 1, 'scale', 'auto'],
        'kernel': ['rbf']
    }

    random_search = RandomizedSearchCV(
        svm, param_distributions=param_dist_random, n_iter=50,
        cv=KFCV, scoring="f1", n_jobs=-1, random_state=42
    )

    start = time_module.time()
    random_search.fit(X_train, y_train)
    elapsed_random = time_module.time() - start

    random_search_output = os.path.join(output_dir, f"RandomSearchCV_results_{config_name}.tsv")
    pd.DataFrame(random_search.cv_results_).to_csv(random_search_output, sep="\t", index=False)

    best_params_random = random_search.best_params_

    if isinstance(best_params_random['gamma'], (int, float)):
        gamma_grid = [best_params_random['gamma'] * 0.1, best_params_random['gamma'], best_params_random['gamma'] * 10]
    else:
        gamma_grid = [best_params_random['gamma']]

    param_grid = {
        'C': [best_params_random['C'] * 0.1, best_params_random['C'], best_params_random['C'] * 2],
        'gamma': gamma_grid,
        'kernel': ['rbf']
    }

    grid_search = GridSearchCV(svm, param_grid=param_grid, cv=KFCV, scoring="f1", n_jobs=-1)
    start = time_module.time()
    grid_search.fit(X_train, y_train)
    elapsed_grid = time_module.time() - start

    grid_search_output = os.path.join(output_dir, f"GridSearchCV_results_{config_name}.tsv")
    pd.DataFrame(grid_search.cv_results_).to_csv(grid_search_output, sep="\t", index=False)

    best_model = grid_search.best_estimator_
    y_pred = best_model.predict(X_test)
    y_proba = best_model.decision_function(X_test)

    f1 = f1_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_proba)

    return {
        "Configuration": config_name,
        "Best_C": grid_search.best_params_['C'],
        "Best_gamma": grid_search.best_params_['gamma'],
        "F1": f1,
        "ROC_AUC": auc,
        "RandomSearch_time_sec": round(elapsed_random, 2),
        "GridSearch_time_sec": round(elapsed_grid, 2)
    }

# ==============================================================================
#  Section 5 — Run Tuning + Evaluation for Each Configuration
# ==============================================================================

results = []

for i, (config_name, transformed_df) in enumerate(configurations.items(), 1):
    X = transformed_df.loc[transformed_df.index.intersection(pfi_species_list)].T

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    result = tune_and_evaluate_svm_rbf(X_train, X_test, y_train, y_test, config_name, OUTPUT_DIR)
    results.append(result)

# ==============================================================================
#  Section 6 — Save Summary Results Table
# ==============================================================================

results_df = pd.DataFrame(results)

summary_output = os.path.join(OUTPUT_DIR, "Summary_Table_all_configurations.csv")
results_df.to_csv(summary_output, index=False)

# ==============================================================================
#  Section 7 — Generate Figures (One Per Technique)
#
# The summary table just saved is reopened from disk here, and a separate F1/ROC-AUC figure is generated and saved for each configuration.
# ==============================================================================

summary_output = os.path.join(OUTPUT_DIR, "Summary_Table_all_configurations.csv")
plot_df = pd.read_csv(summary_output)

for _, row in plot_df.iterrows():
    config_name = row["Configuration"]
    f1_val = row["F1"]
    auc_val = row["ROC_AUC"]

    fig, ax = plt.subplots(figsize=(4.5, 5), dpi=100)

    bar_labels = ["F1", "ROC-AUC"]
    bar_values = [f1_val, auc_val]
    bar_colors = ["#3182bd", "#e6550d"]

    bars = ax.bar(bar_labels, bar_values, color=bar_colors, alpha=0.85,
                   edgecolor="black", linewidth=0.8, width=0.55)

    for bar, val in zip(bars, bar_values):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.01, f"{val:.3f}",
                ha="center", va="bottom", fontsize=11, fontweight="bold")

    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score", fontweight="bold")
    ax.set_title(f"{config_name}", fontweight="bold", fontsize=11)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    plt.tight_layout()

    fig_png = os.path.join(OUTPUT_DIR, f"Figure_{config_name}.png")
    fig_pdf = os.path.join(OUTPUT_DIR, f"Figure_{config_name}.pdf")
    fig_tiff = os.path.join(OUTPUT_DIR, f"Figure_{config_name}.tiff")

    fig.savefig(fig_png, dpi=600, bbox_inches="tight")
    fig.savefig(fig_pdf, bbox_inches="tight")
    fig.savefig(fig_tiff, dpi=600, bbox_inches="tight")

    plt.close(fig)

# ==============================================================================
#  Section 8 — Combined Overview Figure (All Configurations Side-by-Side)
# ==============================================================================

summary_output = os.path.join(OUTPUT_DIR, "Summary_Table_all_configurations.csv")
combined_plot_df = pd.read_csv(summary_output)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), dpi=100)

metrics = ["F1", "ROC_AUC"]
metric_labels = ["Test Set F1 Score", "Test Set ROC-AUC"]
bar_colors = ["#3182bd", "#e6550d", "#31a354", "#756bb1", "#636363"]

x_labels = combined_plot_df["Configuration"].tolist()
x_pos = np.arange(len(x_labels))

for ax, metric, metric_label in zip(axes, metrics, metric_labels):
    values = combined_plot_df[metric].values
    bars = ax.bar(x_pos, values, color=bar_colors, alpha=0.85, edgecolor="black", linewidth=0.8, width=0.6)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.01, f"{val:.3f}",
                ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel(metric_label, fontweight="bold")
    ax.set_title(f"{metric_label} - All Configurations", fontweight="bold", fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

plt.tight_layout()

combined_png = os.path.join(OUTPUT_DIR, "Figure_ALL_configurations_combined.png")
combined_pdf = os.path.join(OUTPUT_DIR, "Figure_ALL_configurations_combined.pdf")
combined_tiff = os.path.join(OUTPUT_DIR, "Figure_ALL_configurations_combined.tiff")

fig.savefig(combined_png, dpi=600, bbox_inches="tight")
fig.savefig(combined_pdf, bbox_inches="tight")
fig.savefig(combined_tiff, dpi=600, bbox_inches="tight")

plt.close(fig)
