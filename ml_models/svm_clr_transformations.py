#! /lustrehome/babluuniba2022/miniconda3/envs/ml/bin/python

# ==============================================================================
#
# This script evaluates whether normalization strategy (Raw counts, TSS, CLR
# additive-pseudocount, and multiplicative zero-replacement) and CLR
# pseudocount value (1.0, 0.5, 0.1) influence the predictive performance of
# our best-performing model (SVM-RBF, PFI-selected features), using the
# Six configurations are evaluated:
# 1. Raw counts (no normalization)
# 2. TSS (Total Sum Scaling / relative abundance)
# 3. CLR, pseudocount = 1.0 (value used in the primary manuscript analysis)
# 4. CLR, pseudocount = 0.5
# 5. CLR, pseudocount = 0.1
# 6. CLR with multiplicative zero-replacement (scikit-bio; zCompositions-style)
#
# All input data is loaded once at the start. All outputs (summary table,
# one figure per technique, and a combined overview figure) are saved to a
# single, clearly named output folder.
# ==============================================================================

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from skbio.stats.composition import multiplicative_replacement, clr
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import f1_score, roc_auc_score

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.size'] = 11

OUTPUT_DIR = "rerun_normalizations"
os.makedirs(OUTPUT_DIR, exist_ok=True)

FIXED_PARAMS = {"C": 20, "gamma": 0.0001, "kernel": "rbf"}


meta_data = pd.read_csv("otu_data/metadata_for_7452.tsv", sep="\t")

training_corrected = pd.read_csv("otu_data/corrected_training_7452.csv", sep="\t", header=0, index_col=0, engine="c")

pfi_dir = "/lustrehome/babluuniba2022/machine_learning/prev_model/model_dev"
pfi_file = os.path.join(pfi_dir, "PERMUTATION_IMPORTANCE_BACTERIAL_bacterial_20pct.tsv")
pfi_results = pd.read_csv(pfi_file, sep="\t")

SPECIES_COLUMN = pfi_results.columns[0]
IMPORTANCE_COLUMN = pfi_results.columns[1]

pfi_species_list = pfi_results[pfi_results[IMPORTANCE_COLUMN] > 0][SPECIES_COLUMN].tolist()


meta_data_indexed = meta_data.set_index("Run_ID")
sample_ids = training_corrected.columns.tolist()
health_labels = meta_data_indexed.loc[sample_ids, "Health_status"]

label_map = {"HEALTHY": 1, "NON-HEALTHY": 0}  # label mapping for our target variable
y = health_labels.map(label_map).values

if pd.isna(y).any():
    raise ValueError("Some Health_status values did not match 'HEALTHY' or 'NON-HEALTHY' exactly.")


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


def normalize_clr_multiplicative_replacement(count_df):
    """CONFIGURATION 6: Centered Log-Ratio transformation with multiplicative
    zero-replacement (scikit-bio implementation, zCompositions-style
    Bayesian-multiplicative approach), an alternative zero-handling method
    to additive pseudocounts."""
    counts = count_df.T.astype(float).values
    sample_totals = counts.sum(axis=1, keepdims=True)
    proportions = counts / sample_totals

    replaced = multiplicative_replacement(proportions)
    clr_transformed = clr(replaced)

    clr_df = pd.DataFrame(
        clr_transformed,
        index=count_df.columns,
        columns=count_df.index
    ).T
    return clr_df


configurations = {
    "Raw_Counts": normalize_raw_counts(training_corrected),
    "TSS_Relative_Abundance": normalize_tss(training_corrected),
    "CLR_pseudocount_1.0": normalize_clr(training_corrected, 1.0),
    "CLR_pseudocount_0.5": normalize_clr(training_corrected, 0.5),
    "CLR_pseudocount_0.1": normalize_clr(training_corrected, 0.1),
    "CLR_multiplicative_replacement": normalize_clr_multiplicative_replacement(training_corrected),
}


results = []

for i, (config_name, transformed_df) in enumerate(configurations.items(), 1):
    X = transformed_df.loc[transformed_df.index.intersection(pfi_species_list)].T

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    print(f"Training set size: {X_train.shape[0]} samples, Test set size: {X_test.shape[0]} samples")
    print(f"Evaluating configuration {i}/{len(configurations)}: {config_name}")
    model = SVC(probability=True, random_state=42, **FIXED_PARAMS)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)[:, 1]

    f1 = f1_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_proba)

    results.append({
        "Configuration": config_name,
        "C": FIXED_PARAMS["C"],
        "gamma": FIXED_PARAMS["gamma"],
        "F1": f1,
        "ROC_AUC": auc,
        "N_features": X.shape[1]
    })
    print(f"Configuration: {config_name}, F1: {f1:.4f}, ROC-AUC: {auc:.4f}, Features used: {X.shape[1]}")



results_df = pd.DataFrame(results)

summary_output = os.path.join(OUTPUT_DIR, "rerun_normalizations_summary.csv")
results_df.to_csv(summary_output, index=False)
summary_output = os.path.join(OUTPUT_DIR, "rerun_normalizations_summary.csv")
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

    fig_png = os.path.join(OUTPUT_DIR, f"rerun_norm_{config_name}.png")
    fig_pdf = os.path.join(OUTPUT_DIR, f"rerun_norm_{config_name}.pdf")
    fig_tiff = os.path.join(OUTPUT_DIR, f"rerun_norm_{config_name}.tiff")

    fig.savefig(fig_png, dpi=600, bbox_inches="tight")
    fig.savefig(fig_pdf, bbox_inches="tight")
    fig.savefig(fig_tiff, dpi=600, bbox_inches="tight")

    plt.close(fig)

summary_output = os.path.join(OUTPUT_DIR, "rerun_normalizations_summary.csv")
combined_plot_df = pd.read_csv(summary_output)

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=100)

metrics = ["F1", "ROC_AUC"]
metric_labels = ["Test Set F1 Score", "Test Set ROC-AUC"]
bar_colors = ["#3182bd", "#e6550d", "#31a354", "#756bb1", "#636363", "#c994c7"]

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

combined_png = os.path.join(OUTPUT_DIR, "rerun_normalizations_ALL_combined.png")
combined_pdf = os.path.join(OUTPUT_DIR, "rerun_normalizations_ALL_combined.pdf")
combined_tiff = os.path.join(OUTPUT_DIR, "rerun_normalizations_ALL_combined.tiff")

fig.savefig(combined_png, dpi=600, bbox_inches="tight")
fig.savefig(combined_pdf, bbox_inches="tight")
fig.savefig(combined_tiff, dpi=600, bbox_inches="tight")

plt.close(fig)

print(f"All outputs saved to: {os.path.abspath(OUTPUT_DIR)}")
