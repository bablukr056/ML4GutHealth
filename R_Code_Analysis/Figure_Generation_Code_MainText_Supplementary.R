# ==============================================================================
# Publication Figure Generation Script
# A Data-Driven Universal Gut Microbiome Health Assessment: 
#
# This script generates all Main Text and Supplementary figures reported in
# the manuscript from the processed model-performance, ROC, feature-importance,
# and metadata tables. Figures are organized and numbered exactly as cited in
# the text (Main Text Figures 1-5; Supplementary Figures 1 and 6).
#
# Output:
#   PUBLICATION_2026/Main_Figures/          -> Figures 1-5
#   PUBLICATION_2026/Supplementary_Figures/  -> Supplementary Figures 1, 6
# ==============================================================================

# ---- Working directory ------------------------------------------------------
setwd("/Users/bablu/OneDrive - Università degli Studi di Milano/CODE/Bablu_ML_Result/")

# ---- Libraries ----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(pheatmap)
library(pracma)
library(scales)
library(grid)

# ---- Output directories --------------------------------------------------------
dir_main <- "PUBLICATION_2026/Main_Figures"
dir_supp <- "PUBLICATION_2026/Supplementary_Figures"
dir.create(dir_main, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_supp, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# SECTION 1 : MAIN TEXT FIGURES
# ==============================================================================

# ------------------------------------------------------------------------------
# Figure 1B - Training Dataset: Phenotype Distribution Across Studies
# ------------------------------------------------------------------------------
bioproject_summary <- read.csv("training_metadata_100K_feb_2026.tsv", sep = "\t")

phenotype_summary <- bioproject_summary %>%
  group_by(Phenotype, Full_Name) %>%
  summarise(Sample_Count = n(), BioProjects = n_distinct(BioProject), .groups = "drop") %>%
  mutate(
    Full_Name = case_when(
      Full_Name == "Atherosclerotic cardiovascular disease" ~ "Atherosclerotic\ncardiovascular disease",
      Full_Name == "Age-related macular degeneration"       ~ "Age-related macular\ndegeneration",
      Full_Name == "Gestational diabetes mellitus"          ~ "Gestational diabetes\nmellitus",
      TRUE ~ Full_Name
    ),
    Phenotype = case_when(Phenotype == "Obesity" ~ "OB", TRUE ~ Phenotype),
    Phenotype_Label = paste0(Full_Name, " (", Phenotype, ")")
  )

plot_data <- phenotype_summary %>%
  select(-Phenotype) %>%
  mutate(
    Phenotype_Label = ifelse(Phenotype_Label == "Healthy (Healthy)", "Healthy (HC)", Phenotype_Label),
    Percentage = round(Sample_Count / sum(Sample_Count) * 100, 1),
    Phenotype_Label_Display = Phenotype_Label
  )

pheno_colors <- c(
  "Healthy (HC)" = "#4ECDC4", "Colorectal Cancer (CRC)" = "#E63946",
  "Antibiotic exposer (ABx)" = "#1D3557", "Parkinson (PD)" = "#9B59B6",
  "Obesity (OB)" = "#2ECC71", "Crohn Disease (CD)" = "#F39C12",
  "Ulcerative Colitis (UC)" = "#E67E22",
  "Atherosclerotic\ncardiovascular disease (ACVD)" = "#95A5A6",
  "Acute Pancreatitis (AP)" = "#16A085", "Liver cirrhosis (LV)" = "#8E44AD",
  "Ankylosing spondylitis (AS)" = "#27AE60", "Gastric Cancer (GC)" = "#E91E63",
  "Age-related macular\ndegeneration (AMD)" = "#D68910",
  "Gestational diabetes\nmellitus (GDM)" = "#A04000"
)

plot_data <- plot_data %>%
  mutate(
    Label_With_Info = sprintf("%s\nn=%d", Phenotype_Label, Sample_Count),
    Label_With_Info = str_wrap(Label_With_Info, width = 35),
    Label_With_Info = fct_reorder(Label_With_Info, Sample_Count, .desc = FALSE)
  )

max_sample_count <- max(plot_data$Sample_Count)
x_limit <- max_sample_count * 1.28

fig1B <- ggplot(plot_data, aes(x = Sample_Count, y = Label_With_Info, fill = Phenotype_Label_Display)) +
  geom_bar(stat = "identity", color = "#000000", linewidth = 0.4) +
  geom_text(
    aes(x = Sample_Count, y = Label_With_Info, label = sprintf("%d \u2192 %.1f%%", BioProjects, Percentage)),
    hjust = -0.1, vjust = 0.5, color = "#000000", fontface = "bold", size = 6, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = pheno_colors) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, x_limit), breaks = seq(0, max_sample_count, by = 500)) +
  labs(x = "# of Samples", y = "",
       title = "7,452 Stool-Derived Shotgun Metagenomes Across 13 Phenotypes\nfrom 32 Independent Studies") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "#1a1a1a", margin = margin(b = 4), lineheight = 1.3),
    axis.text.y = element_text(size = 14, color = "#000000", face = "bold", hjust = 1, lineheight = 1.1),
    axis.text.x = element_text(size = 14, color = "#000000", face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 10)),
    axis.line.x = element_line(color = "#000000", linewidth = 0.6),
    axis.line.y = element_blank(),
    axis.ticks.x = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(0.15, "cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 20, 12, 15),
    legend.position = "none"
  )

print(fig1B)
ggsave(file.path(dir_main, "Fig1B_training_dataset_phenotype_distribution.png"), fig1B, width = 14, height = 9, dpi = 600, bg = "white")
ggsave(file.path(dir_main, "Fig1B_training_dataset_phenotype_distribution.pdf"), fig1B, width = 14, height = 9, device = "pdf", bg = "white")


# ------------------------------------------------------------------------------
# Figure 2A - F1 Score: Feature Selection Methods Performance (Test Dataset)
# ------------------------------------------------------------------------------
metrics_file <- read.csv("perfomance_on_test_dataset_combined_for_all_models_prev_0_to_25.tsv", sep = "\t")

metrics_file <- metrics_file %>%
  filter(Prevalence == 20, grepl("^(SVM_RBF|SVM_LIN|RF|LogReg)", Model)) %>%
  mutate(
    BaseModel = case_when(
      grepl("^SVM_RBF", Model) ~ "SVM-RBF",
      grepl("^SVM_LIN", Model) ~ "SVM-Linear",
      grepl("^RF", Model)      ~ "RF",
      grepl("^LogReg", Model)  ~ "LogReg-ElasticNet"
    ),
    FeatureSelection = case_when(
      Model %in% c("SVM_RBF", "SVM_LIN", "RF", "LogReg") ~ paste0("All Features\n(n=", Features, ")"),
      grepl("_LASSO$", Model) ~ paste0("LASSO\n(n=", Features, ")"),
      grepl("_RFECV$", Model) ~ paste0("RFECV\n(n=", Features, ")"),
      grepl("_PFI$", Model)   ~ paste0("PFI\n(n=", Features, ")"),
      TRUE ~ paste0("All Features\n(n=", Features, ")")
    ),
    FeatureSelection = factor(FeatureSelection, levels = unique(FeatureSelection[order(Features, decreasing = TRUE)])),
    F1_Percent = F1.Score * 100
  )

custom_colors <- c("SVM-RBF" = "#E63946", "SVM-Linear" = "#1D3557", "RF" = "#2A9D8F", "LogReg-ElasticNet" = "#F4A261")

fig2A <- ggplot(metrics_file, aes(x = FeatureSelection, y = F1_Percent, fill = BaseModel)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75, color = "#000000", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f%%", F1_Percent), group = BaseModel),
            position = position_dodge(width = 0.8), vjust = -0.6, size = 3.8, color = "#000000", fontface = "bold") +
  labs(title = "Model Performance (F1 Score) on Test Dataset: n = 1,491 (798 non-healthy, 693 healthy)",
       x = "Feature Selection Method (Number of Features)", y = "F1 Score (%)", fill = "Model") +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), breaks = seq(0, 100, by = 10), limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0, color = "#1a1a1a", margin = margin(b = 3)),
    axis.title.x = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    axis.text.x = element_text(size = 12, face = "bold", color = "#000000", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    legend.text = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 12)),
    legend.key.size = unit(0.5, "cm"),
    legend.key.width = unit(1.3, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 10, b = 0),
    plot.margin = margin(12, 12, 10, 12),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(fig2A)
ggsave(file.path(dir_main, "Fig2A_F1_score_feature_selection_comparison.png"), fig2A, width = 11, height = 6, dpi = 600, bg = "white")
ggsave(file.path(dir_main, "Fig2A_F1_score_feature_selection_comparison.pdf"), fig2A, width = 11, height = 6, device = "pdf", bg = "white")


# ------------------------------------------------------------------------------
# Figure 2B - ROC Curves Across Feature Selection Methods (Test Dataset)
# ------------------------------------------------------------------------------
roc_data <- read.csv("roc_data_test_dataset_20pct.csv", stringsAsFactors = FALSE)

family_colors <- c("SVM_RBF" = "#E63946", "SVM_LIN" = "#1D3557", "RF" = "#2A9D8F", "LogReg" = "#F4A261")

roc_data <- roc_data %>%
  mutate(
    BaseModel = case_when(
      grepl("^SVM_RBF", Model) ~ "SVM_RBF",
      grepl("^SVM_LIN", Model) ~ "SVM_LIN",
      grepl("^RF", Model)      ~ "RF",
      grepl("^LogReg", Model)  ~ "LogReg",
      TRUE ~ Model
    ),
    Feature_Selection = case_when(
      Model %in% c("SVM_RBF", "SVM_LIN", "RF", "LogReg") ~ "All Features",
      grepl("_LASSO$", Model) ~ "LASSO",
      grepl("_RFECV$", Model) ~ "RFECV",
      grepl("_PFI$", Model)   ~ "PFI",
      TRUE ~ Feature_Selection
    ),
    Feature_Selection = factor(Feature_Selection, levels = c("All Features", "LASSO", "PFI", "RFECV")),
    Feature_Selection_Label = case_when(
      Feature_Selection == "All Features" ~ "All Features (n=4,022)",
      Feature_Selection == "LASSO"        ~ "LASSO (n=828)",
      Feature_Selection == "PFI"          ~ "PFI (n=3,127)",
      Feature_Selection == "RFECV"        ~ "RFECV (n=3,472)"
    ),
    Feature_Selection_Label = factor(Feature_Selection_Label,
                                     levels = c("All Features (n=4,022)", "LASSO (n=828)", "PFI (n=3,127)", "RFECV (n=3,472)")),
    DisplayModel = case_when(
      grepl("^SVM_RBF", Model) ~ "SVM-RBF",
      grepl("^SVM_LIN", Model) ~ "SVM-Linear",
      grepl("^RF", Model)      ~ "RF",
      grepl("^LogReg", Model)  ~ "LogReg (ElasticNet)"
    )
  )

auc_data <- roc_data %>%
  group_by(Model, DisplayModel, Feature_Selection, Feature_Selection_Label, BaseModel) %>%
  arrange(FPR) %>%
  summarise(AUC = trapz(FPR, TPR), .groups = "drop") %>%
  group_by(Feature_Selection_Label) %>%
  arrange(desc(AUC)) %>%
  mutate(rank = row_number(), AUC_percent = AUC * 100,
         auc_label = sprintf("%s (%.1f%%)", DisplayModel, AUC_percent),
         label_x = 0.39, label_y = 0.38 - (rank - 1) * 0.085) %>%
  ungroup()

box_data <- auc_data %>% group_by(Feature_Selection_Label) %>% slice(1) %>% ungroup()

fig2B <- ggplot(roc_data, aes(x = FPR, y = TPR, color = BaseModel)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.55) +
  geom_line(linewidth = 1.4, alpha = 0.92) +
  facet_grid(. ~ Feature_Selection_Label) +
  geom_rect(data = box_data, aes(xmin = 0.36, xmax = 1.00, ymin = 0.00, ymax = 0.41),
            inherit.aes = FALSE, fill = "white", color = "#bbbbbb", linewidth = 0.5) +
  geom_text(data = auc_data, aes(x = label_x, y = label_y, label = auc_label, color = BaseModel),
            size = 4, fontface = "bold", hjust = 0, show.legend = FALSE) +
  labs(title = "ROC Curve Comparison of All ML Models on Test Dataset (n = 1,491)",
       x = "False Positive Rate", y = "True Positive Rate", color = "Model") +
  scale_color_manual(values = family_colors,
                     labels = c("SVM_RBF" = "SVM-RBF", "SVM_LIN" = "SVM-Linear", "RF" = "RF", "LogReg" = "LogReg (ElasticNet)")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01), labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01), labels = function(x) sprintf("%.1f", x)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "#1a1a1a", margin = margin(b = 8), lineheight = 1.3),
    axis.title.x = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 8)),
    axis.text = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.12, "cm"),
    strip.text = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 5, b = 5)),
    strip.background = element_rect(fill = "white", color = "#000000", linewidth = 0.7),
    panel.spacing = unit(0.8, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 8)),
    legend.text = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    legend.key.width = unit(1.4, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 6, b = 4),
    plot.margin = margin(12, 12, 10, 12)
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 2.5)))

print(fig2B)
ggsave(file.path(dir_main, "Fig2B_roc_curves_test_dataset.png"), fig2B, width = 14, height = 6, dpi = 600, bg = "white")
ggsave(file.path(dir_main, "Fig2B_roc_curves_test_dataset.pdf"), fig2B, width = 18, height = 6, device = "pdf", bg = "white")


# ------------------------------------------------------------------------------
# Figure 2C - Class-Wise Model Performance Across Feature Selection Methods (Test)
# ------------------------------------------------------------------------------
combined_df <- read.csv("classification_report_for_all_feature_selection_20pct_models.csv", sep = "\t")
combined_df <- combined_df %>% mutate(Prevalence = as.integer(gsub("pct", "", Prevalence)))

feature_counts <- c("All_Features" = 4022, "LASSO" = 828, "PFI" = 3127, "RFE" = 3472)

model_families <- list(
  SVM_RBF = c("SVM_RBF", "SVM_RBF_PFI", "SVM_RBF_RFECV", "SVM_RBF_LASSO"),
  SVM_LIN = c("SVM_LIN", "SVM_LIN_PFI", "SVM_LIN_RFECV", "SVM_LIN_LASSO"),
  RF      = c("RF", "RF_LASSO", "RF_PFI", "RF_RFECV"),
  LogReg  = c("LogReg", "LogReg_PFI", "LogReg_RFECV", "LogReg_LASSO")
)

model_to_family <- unlist(lapply(names(model_families), function(fam) {
  setNames(rep(fam, length(model_families[[fam]])), model_families[[fam]])
}))

df_combined <- combined_df %>%
  filter(Prevalence == 20, Dataset == "test", Class_Label %in% c("0", "1")) %>%
  mutate(
    Class_Label_Desc = case_when(Class_Label == "0" ~ "Non-Healthy", Class_Label == "1" ~ "Healthy", TRUE ~ Class_Label),
    Model_Display = gsub("_RFECV", "_RFE", Model),
    Model_Display = gsub("SVM_RBF", "SVM-RBF", Model_Display),
    Model_Display = gsub("SVM_LIN", "SVM-Linear", Model_Display),
    Model_Display = gsub("_", " ", Model_Display),
    Model_Family = model_to_family[Model],
    F1_Percent = F1_Score * 100,
    Feature_Selection_Label = case_when(
      Feature_Selection == "All_Features" ~ paste0("All Features (n=", feature_counts["All_Features"], ")"),
      Feature_Selection == "LASSO"        ~ paste0("LASSO (n=", feature_counts["LASSO"], ")"),
      Feature_Selection == "PFI"          ~ paste0("PFI (n=", feature_counts["PFI"], ")"),
      Feature_Selection == "RFE"          ~ paste0("RFECV (n=", feature_counts["RFE"], ")")
    ),
    Feature_Selection_Label = factor(Feature_Selection_Label, levels = c(
      paste0("All Features (n=", feature_counts["All_Features"], ")"),
      paste0("LASSO (n=", feature_counts["LASSO"], ")"),
      paste0("PFI (n=", feature_counts["PFI"], ")"),
      paste0("RFECV (n=", feature_counts["RFE"], ")")
    ))
  ) %>%
  group_by(Feature_Selection_Label) %>%
  mutate(Model_Display = factor(Model_Display, levels = unique(Model_Display[order(Model_Family, Model_Display)]))) %>%
  ungroup()

fig2C <- ggplot(df_combined, aes(x = Model_Display, y = F1_Percent, fill = Class_Label_Desc)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.7, color = "#000000", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", F1_Percent)), position = position_dodge(width = 0.75),
            vjust = -0.7, size = 4.4, color = "#000000", fontface = "bold") +
  facet_wrap(~ Feature_Selection_Label, nrow = 2, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = c("Non-Healthy" = "#E69F00", "Healthy" = "#56B4E9"),
                    labels = c("Non-Healthy" = "Non-Healthy (n=798)", "Healthy" = "Healthy (n=693)")) +
  labs(title = "Class-Wise Model Performance (F1 Score) on Test Dataset (n=1,491)", x = "", y = "F1 Score (%)", fill = "Class") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20)), limits = c(0, 100), breaks = seq(0, 100, by = 20), labels = function(x) paste0(x, "%")) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0, color = "#1a1a1a", margin = margin(b = 4), lineheight = 1.25),
    axis.title.y = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    axis.text.x = element_text(size = 11, face = "bold", color = "#000000", angle = 0, hjust = 0.5, vjust = 1),
    axis.text.y = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks = element_line(color = "#000000", linewidth = 0.5),
    strip.text = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 5, b = 5)),
    strip.background = element_rect(fill = "white", color = "#000000", linewidth = 0.7),
    panel.spacing = unit(1.1, "lines"),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "top",
    legend.title = element_text(size = 14, face = "bold", color = "#000000", margin = margin(r = 10)),
    legend.text = element_text(size = 14, face = "bold", color = "#000000", margin = margin(r = 12)),
    legend.key.width = unit(2, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 5, b = 10),
    plot.margin = margin(14, 14, 12, 14)
  )

print(fig2C)
ggsave(file.path(dir_main, "Fig2C_classwise_F1_test_dataset_all_feature_selections.png"), fig2C, width = 14, height = 8, dpi = 600, bg = "white")
ggsave(file.path(dir_main, "Fig2C_classwise_F1_test_dataset_all_feature_selections.pdf"), fig2C, width = 14, height = 8, device = "pdf", bg = "white")


# ------------------------------------------------------------------------------
# Figure 3A - Top 50 Species by Permutation Feature Importance (PFI)
# ------------------------------------------------------------------------------
pfi <- read_tsv("PERMUTATION_IMPORTANCE_BACTERIAL_bacterial_20pct.tsv")
pfi <- pfi %>% filter(Importance_Mean > 0)
pfi_top50 <- pfi %>% arrange(desc(Importance_Mean)) %>% slice_head(n = 50)

fig3A <- ggplot(pfi_top50, aes(x = Importance_Mean, y = reorder(Feature, Importance_Mean), fill = Importance_Mean)) +
  geom_col(color = "#000000", width = 0.75, linewidth = 0.35) +
  geom_errorbar(aes(xmin = pmax(0, Importance_Mean - Importance_Std), xmax = Importance_Mean + Importance_Std),
                width = 0.4, linewidth = 0.5, color = "#000000") +
  geom_text(aes(label = sprintf("%.4f", Importance_Mean)), hjust = -0.15, size = 3.2, color = "#000000", fontface = "bold") +
  scale_fill_gradientn(colors = c("#3182bd", "#FFFFFF", "#e34a33"), name = "Importance Score") +
  labs(title = "Top 50 Bacterial Species by SVM-RBF-PFI", x = "Mean Permutation Importance Score (\u00b1 Std)", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "#1a1a1a", margin = margin(b = 6), lineheight = 1.3),
    axis.title.x = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 10)),
    axis.text.y = element_text(size = 12, face = "bold.italic", color = "#000000"),
    axis.text.x = element_text(size = 11, face = "bold", color = "#000000"),
    axis.line = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks = element_line(color = "#000000", linewidth = 0.4),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major.x = element_line(color = "white", linewidth = 0.25),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(14, 14, 12, 14),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 11, face = "bold", color = "#000000", margin = margin(r = 8)),
    legend.text = element_text(size = 10, face = "bold", color = "#000000"),
    legend.key.width = unit(3.0, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 5, b = 10)
  )

print(fig3A)
ggsave(file.path(dir_main, "Fig3A_permutation_importance_top50_svm_rbf_pfi.png"), fig3A, width = 12, height = 14, dpi = 600, bg = "white")
ggsave(file.path(dir_main, "Fig3A_permutation_importance_top50_svm_rbf_pfi.pdf"), fig3A, width = 12, height = 14, device = "pdf", bg = "white")


# ------------------------------------------------------------------------------
# Figure 3B - MaAsLin2: Significant Associations of Top 50 PFI Features
# ------------------------------------------------------------------------------
sig_df <- read.delim("top_50_pfi_maaslin_significant.tsv", sep = "\t")

maaslin2_heatmap <- function(df, title = "Top 50 Microbial Features: Significant Associations\n(-log(qval) * sign(coef))",
                             cell_value = "qval", border_color = "#000000",
                             color = colorRampPalette(c("#3182bd", "#FFFFFF", "#e34a33")), first_n = 50) {
  
  if (!is.na(first_n) & first_n > 0 & first_n < nrow(df)) {
    df <- if (cell_value == "coef") df[order(-abs(df[[cell_value]])), ] else df[order(df[[cell_value]]), ]
    df_sub <- df[1:first_n, ]
    for (idx in seq(first_n, nrow(df))) {
      if (length(unique(df_sub$feature)) == first_n) break
      df_sub <- df[1:idx, ]
    }
    df <- df[df$feature %in% df_sub$feature, ]
  }
  if (nrow(df) < 2) { message("No associations to plot."); return(NULL) }
  
  metadata <- df$metadata; data <- df$feature; dfvalue <- df$value
  value <- if (cell_value == "pval") pmax(-20, pmin(20, -log(df$pval) * sign(df$coef)))
  else if (cell_value == "qval") pmax(-20, pmin(20, -log(df$qval) * sign(df$coef)))
  else df$coef
  
  verbose_metadata <- c(); metadata_multi_level <- c()
  for (i in unique(metadata)) {
    levels <- unique(df$value[df$metadata == i])
    if (length(levels) > 1) {
      metadata_multi_level <- c(metadata_multi_level, i)
      verbose_metadata <- c(verbose_metadata, levels)
    } else {
      verbose_metadata <- c(verbose_metadata, if (i == "Health_status") "Healthy" else i)
    }
  }
  
  n <- length(unique(data)); m <- length(unique(verbose_metadata))
  if (n < 2 || m < 2) { message("Not enough features/metadata for heatmap."); return(NULL) }
  
  a <- as.data.frame(matrix(0, nrow = n, ncol = m))
  rownames(a) <- unique(data); colnames(a) <- unique(verbose_metadata)
  
  for (i in seq_len(nrow(df))) {
    cm <- metadata[i]
    if (cm %in% metadata_multi_level) cm <- dfvalue[i] else if (cm == "Health_status") cm <- "Healthy"
    if (abs(a[as.character(data[i]), as.character(cm)]) > abs(value[i])) next
    a[as.character(data[i]), as.character(cm)] <- value[i]
  }
  
  max_value <- max(abs(a), na.rm = TRUE); if (max_value == 0) max_value <- 1
  breaks <- seq(-max_value, max_value, length.out = 101)
  
  p <- pheatmap::pheatmap(
    a, cellwidth = 20, cellheight = 12, main = title, fontsize = 13, fontsize_row = 9, fontsize_col = 11,
    border = TRUE, border_color = border_color, show_rownames = TRUE, show_colnames = TRUE, scale = "none",
    cluster_rows = TRUE, cluster_cols = TRUE,
    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
    legend = TRUE,
    legend_breaks = seq(-ceiling(max_value), ceiling(max_value), by = max(1, ceiling(max_value / 5))),
    legend_labels = sprintf("%.1f", seq(-ceiling(max_value), ceiling(max_value), by = max(1, ceiling(max_value / 5)))),
    color = color(length(breaks)), breaks = breaks, treeheight_row = 50, treeheight_col = 50,
    display_numbers = matrix(ifelse(a > 0, "+", ifelse(a < 0, "-", "")), nrow(a)),
    number_color = "#000000", fontsize_number = 8, angle_col = "90", silent = FALSE
  )
  message("Heatmap value range: ", round(min(a, na.rm = TRUE), 2), " to ", round(max(a, na.rm = TRUE), 2))
  return(p)
}

fig3B <- maaslin2_heatmap(df = sig_df, cell_value = "qval", first_n = 50)

png(file.path(dir_main, "Fig3B_maaslin2_heatmap_top50_features.png"), width = 4000, height = 4500, res = 400, type = "cairo")
grid.newpage()
grid.draw(fig3B$gtable)
dev.off()


# ------------------------------------------------------------------------------
# Figure 4A - Validation Dataset: Phenotype Distribution
# ------------------------------------------------------------------------------
bioproject_summary <- read.csv("fina_validaation_metadata_for_publications.tsv", sep = "\t")

phenotype_summary <- bioproject_summary %>%
  group_by(Phenotype, Full_Name) %>%
  summarise(Sample_Count = n(), BioProjects = n_distinct(BioProject), .groups = "drop") %>%
  mutate(
    Full_Name = case_when(Full_Name == "Clostridium difficile infection" ~ "Clostridium difficile\ninfection", TRUE ~ Full_Name),
    Phenotype = case_when(Phenotype == "Obesity" ~ "OB", TRUE ~ Phenotype),
    Phenotype_Label = paste0(Full_Name, " (", Phenotype, ")")
  )

plot_data <- phenotype_summary %>%
  select(-Phenotype) %>%
  mutate(
    Phenotype_Label = ifelse(Phenotype_Label == "Healthy (Healthy)", "Healthy (HC)", Phenotype_Label),
    Percentage = round(Sample_Count / sum(Sample_Count) * 100, 1),
    Phenotype_Label_Display = Phenotype_Label
  )

pheno_colors <- c(
  "Healthy (HC)" = "#4ECDC4", "Colorectal Cancer (CRC)" = "#E63946", "Parkinson (PD)" = "#9B59B6",
  "Obesity (OB)" = "#2ECC71", "Crohn Disease (CD)" = "#F39C12",
  "Atherosclerotic cardiovascular disease (ACVD)" = "#95A5A6",
  "Clostridium difficile\ninfection (CDI)" = "#16A085", "Type 2 diabetes (T2D)" = "#FF1493"
)

plot_data <- plot_data %>%
  mutate(
    Label_With_Info = sprintf("%s\nn=%d", Phenotype_Label, Sample_Count),
    Label_With_Info = str_wrap(Label_With_Info, width = 35),
    Label_With_Info = fct_reorder(Label_With_Info, Sample_Count, .desc = FALSE)
  )

max_sample_count <- max(plot_data$Sample_Count)
x_limit <- max_sample_count * 1.32
total_samples <- sum(plot_data$Sample_Count)
total_bioprojects <- n_distinct(bioproject_summary$BioProject)

fig4A <- ggplot(plot_data, aes(x = Sample_Count, y = Label_With_Info, fill = Phenotype_Label_Display)) +
  geom_bar(stat = "identity", color = "#000000", linewidth = 0.45) +
  geom_text(aes(x = Sample_Count, y = Label_With_Info, label = sprintf("%d \u2192 %.1f%%", BioProjects, Percentage)),
            hjust = -0.08, vjust = 0.5, color = "#000000", fontface = "bold", size = 6, inherit.aes = FALSE) +
  scale_fill_manual(values = pheno_colors) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, x_limit), breaks = seq(0, max_sample_count, by = 50)) +
  labs(x = "# of Samples", y = "",
       title = sprintf("%d Stool-Derived Shotgun Metagenomes Across 7 Phenotypes\nfrom %d Independent Studies", total_samples, total_bioprojects)) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "#1a1a1a", margin = margin(b = 12), lineheight = 1.35),
    axis.text.y = element_text(size = 14, color = "#000000", face = "bold", hjust = 1, lineheight = 1.15),
    axis.text.x = element_text(size = 14, color = "#000000", face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold", color = "#000000", margin = margin(t = 10)),
    axis.line.x = element_line(color = "#000000", linewidth = 0.6),
    axis.line.y = element_blank(),
    axis.ticks.x = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.y = element_blank(),
    axis.ticks.length.x = unit(0.15, "cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#eeeeee", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.8),
    plot.margin = margin(18, 24, 14, 16),
    legend.position = "none"
  )

print(fig4A)
ggsave(file.path(dir_main, "Fig4A_validation_dataset_phenotype_distribution.png"), fig4A, width = 15, height = 9, dpi = 600, bg = "white")
ggsave(file.path(dir_main, "Fig4A_validation_dataset_phenotype_distribution.pdf"), fig4A, width = 15, height = 9, device = "pdf", bg = "white")


# ------------------------------------------------------------------------------
# Figure 4B - SVM-RBF F1 Score Across Feature Selection Methods (Validation)
# ------------------------------------------------------------------------------
metrics_file <- read.csv("Performance_on_validation_dataset_20pct_2026.tsv", sep = "\t")

metrics_file <- metrics_file %>%
  filter(grepl("^SVM_RBF", Model)) %>%
  mutate(
    FeatureSelection = case_when(
      Feature_Selection == "All_Features" ~ "All Features",
      Feature_Selection == "RFE"          ~ "RFECV",
      Feature_Selection == "PFI"          ~ "PFI",
      Feature_Selection == "LASSO"        ~ "LASSO",
      TRUE ~ Feature_Selection
    ),
    FeatureLabel = paste0(FeatureSelection, "\n(n=", Features, ")"),
    F1_Percent = F1.Score * 100
  )

custom_colors <- c("All Features" = "#E63946", "PFI" = "#00B7EB", "RFECV" = "#2A9D8F", "LASSO" = "#E9C46A")

fig4B <- ggplot(metrics_file, aes(x = reorder(FeatureLabel, -F1_Percent), y = F1_Percent, fill = FeatureSelection)) +
  geom_bar(stat = "identity", width = 0.65, color = "#000000", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", F1_Percent)), vjust = -0.7, size = 5, color = "#000000", fontface = "bold") +
  labs(title = "SVM-RBF F1 Score on Validation Dataset: n = 642 (426 non-healthy, 216 healthy)",
       x = "Feature Selection Method", y = "F1 Score (%)", fill = "Models") +
  scale_fill_manual(values = custom_colors,
                    labels = c("All Features" = "SVM-RBF (All Features)", "RFECV" = "SVM-RBF (RFECV)",
                               "PFI" = "SVM-RBF (PFI)", "LASSO" = "SVM-RBF (LASSO)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), limits = c(0, 100), breaks = seq(0, 100, by = 20), labels = function(x) paste0(x, "%")) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0, color = "#1a1a1a", margin = margin(b = 4)),
    axis.title.y = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    axis.title.x = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 10)),
    axis.text.x = element_text(size = 12, face = "bold", color = "#000000"),
    axis.text.y = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks = element_line(color = "#000000", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    legend.text = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 12)),
    legend.key.width = unit(1.5, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(fig4B)
ggsave(file.path(dir_main, "Fig4B_F1_score_validation_svm_rbf_all_methods.png"), fig4B, width = 12, height = 6, dpi = 600, bg = "white")
ggsave(file.path(dir_main, "Fig4B_F1_score_validation_svm_rbf_all_methods.pdf"), fig4B, width = 12, height = 6, device = "pdf", bg = "white")


# ------------------------------------------------------------------------------
# Figure 4C - Class-Wise SVM-RBF Performance (Validation Dataset)
# ------------------------------------------------------------------------------
combined_df <- read.csv("classification_reports_validation_dataset_20pct_2026.tsv", sep = "\t")
combined_df <- combined_df %>% mutate(Prevalence = as.integer(gsub("pct", "", Prevalence)))

feature_counts <- c("All_Features" = 4022, "LASSO" = 828, "PFI" = 3127, "RFE" = 3472)

df_filtered <- combined_df %>%
  filter(Prevalence == 20, Dataset == "validation", Class_Label %in% c("0", "1"),
         Feature_Selection %in% names(feature_counts),
         Model %in% c("SVM_RBF", "SVM_RBF_PFI", "SVM_RBF_LASSO", "SVM_RBF_RFECV")) %>%
  mutate(
    Class_Label_Desc = recode(Class_Label, "0" = "Non-Healthy", "1" = "Healthy"),
    Feature_Selection_Clean = recode(Feature_Selection, "RFE" = "RFECV", "All_Features" = "All Features", .default = Feature_Selection),
    Feature_Selection_Annotated = paste0(Feature_Selection_Clean, "\n(n=", feature_counts[Feature_Selection], ")"),
    Feature_Selection_Annotated = factor(Feature_Selection_Annotated,
                                         levels = paste0(c("All Features", "RFECV", "PFI", "LASSO"), "\n(n=", feature_counts[c("All_Features", "RFE", "PFI", "LASSO")], ")")),
    F1_Percent = F1_Score * 100
  )

fig4C <- ggplot(df_filtered, aes(x = Feature_Selection_Annotated, y = F1_Percent, fill = Class_Label_Desc)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.7, color = "#000000", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", F1_Percent)), position = position_dodge(width = 0.75),
            vjust = -0.7, size = 4.4, color = "#000000", fontface = "bold") +
  scale_fill_manual(values = c("Non-Healthy" = "#E69F00", "Healthy" = "#56B4E9"),
                    labels = c("Non-Healthy" = "Non-Healthy (n=472)", "Healthy" = "Healthy (n=165)")) +
  labs(x = "", y = "F1 Score (%)", fill = "Class",
       title = "Class-Wise SVM-RBF Model Performance (F1 Score) on Validation Dataset (n=642)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20)), limits = c(0, 100), breaks = seq(0, 100, by = 20), labels = function(x) paste0(x, "%")) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "#1a1a1a", margin = margin(b = 4), lineheight = 1.3),
    axis.title.y = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    axis.title.x = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 10)),
    axis.text.x = element_text(size = 12, face = "bold", color = "#000000"),
    axis.text.y = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks = element_line(color = "#000000", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "top",
    legend.title = element_text(size = 14, face = "bold", color = "#000000", margin = margin(r = 10)),
    legend.text = element_text(size = 14, face = "bold", color = "#000000", margin = margin(r = 12)),
    legend.key.width = unit(1.5, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 5, b = 10),
    plot.margin = margin(14, 14, 12, 14)
  )

print(fig4C)
ggsave(file.path(dir_main, "Fig4C_validation_F1_classwise_svm_rbf.png"), fig4C, width = 10, height = 6, dpi = 600, bg = "white")
ggsave(file.path(dir_main, "Fig4C_validation_F1_classwise_svm_rbf.pdf"), fig4C, width = 10, height = 6, device = "pdf", bg = "white")


# ------------------------------------------------------------------------------
# Figure 4D - ROC Curve: SVM-RBF Across All Feature Selections (Validation)
# ------------------------------------------------------------------------------
roc_data <- read.csv("roc_data_validation_dataset_20pct_2026.tsv", sep = "\t", stringsAsFactors = FALSE)
roc_data$Model <- gsub("LR[-_]ElasticNet", "LogReg", roc_data$Model)

model_colors <- c("SVM_RBF" = "#E63946", "SVM_RBF_RFE" = "#2A9D8F", "SVM_RBF_PFI" = "#00B7EB", "SVM_RBF_LASSO" = "#E9C46A")

roc_data <- roc_data %>%
  filter(grepl("^SVM_RBF", Model)) %>%
  mutate(
    Feature_Selection = case_when(
      grepl("_RFE$", Model) ~ "RFECV", grepl("_PFI$", Model) ~ "PFI",
      grepl("_LASSO$", Model) ~ "LASSO", Model == "SVM_RBF" ~ "All Features", TRUE ~ "Other"
    ),
    Feature_Selection = factor(Feature_Selection, levels = c("All Features", "RFECV", "PFI", "LASSO")),
    DisplayModel = case_when(
      grepl("_RFE$", Model) ~ "SVM-RBF (RFECV)", grepl("_PFI$", Model) ~ "SVM-RBF (PFI)",
      grepl("_LASSO$", Model) ~ "SVM-RBF (LASSO)", Model == "SVM_RBF" ~ "SVM-RBF (All Features)", TRUE ~ Model
    )
  )

auc_data <- roc_data %>%
  group_by(Model, DisplayModel, Feature_Selection) %>%
  arrange(FPR) %>%
  summarise(AUC = trapz(FPR, TPR), .groups = "drop") %>%
  arrange(desc(AUC)) %>%
  mutate(rank = row_number(), AUC_percent = AUC * 100,
         auc_label = sprintf("%s (%.1f%%)", DisplayModel, AUC_percent),
         label_x = 0.5, label_y = 0.38 - (rank - 1) * 0.085)

fig4D <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.55) +
  geom_line(linewidth = 1.4, alpha = 0.92) +
  geom_rect(data = auc_data[1, ], aes(xmin = 0.49, xmax = 1.00, ymin = 0.1, ymax = 0.41),
            inherit.aes = FALSE, fill = "white", color = "#bbbbbb", linewidth = 0.5) +
  geom_text(data = auc_data, aes(x = label_x, y = label_y, label = auc_label, color = Model),
            size = 4, fontface = "bold", hjust = 0, show.legend = FALSE) +
  labs(title = "ROC Curve Comparison of SVM-RBF Models on Validation Dataset (n = 642)",
       x = "False Positive Rate", y = "True Positive Rate", color = "Model") +
  scale_color_manual(values = model_colors,
                     labels = c("SVM_RBF" = "SVM-RBF (All Features)", "SVM_RBF_RFE" = "SVM-RBF (RFECV)",
                                "SVM_RBF_PFI" = "SVM-RBF (PFI)", "SVM_RBF_LASSO" = "SVM-RBF (LASSO)")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01), labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01), labels = function(x) sprintf("%.1f", x)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "#1a1a1a", margin = margin(b = 8), lineheight = 1.3),
    axis.title.x = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 8)),
    axis.text = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 8)),
    legend.text = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    legend.key.width = unit(1.4, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 6, b = 4),
    plot.margin = margin(12, 12, 10, 12)
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 2.5)))

print(fig4D)
ggsave(file.path(dir_main, "Fig4D_roc_validation_svm_rbf.png"), fig4D, width = 12, height = 8, dpi = 600, bg = "white")
ggsave(file.path(dir_main, "Fig4D_roc_validation_svm_rbf.pdf"), fig4D, width = 12, height = 8, device = "pdf", bg = "white")


# ------------------------------------------------------------------------------
# Figure 5B - ROC Curves by BioProject (SVM-RBF, Validation Datasets)
# ------------------------------------------------------------------------------
roc_data <- read.delim("roc_data_bioproject_20pct_2026.tsv", sep = "\t", header = TRUE)

roc_data <- roc_data %>%
  filter(grepl("^SVM_RBF", Model)) %>%
  mutate(
    Phenotype_Counts = str_replace_all(Phenotype_Counts, ": ", ":"),
    ModelDisplay = case_when(
      grepl("_RFE", Model) ~ "SVM-RBF (RFECV)", grepl("_PFI", Model) ~ "SVM-RBF (PFI)",
      grepl("_LASSO", Model) ~ "SVM-RBF (LASSO)", Model == "SVM_RBF" ~ "SVM-RBF (All Features)", TRUE ~ Model
    )
  )

model_colors <- c("SVM_RBF" = "#E63946", "SVM_RBF_RFECV" = "#2A9D8F", "SVM_RBF_PFI" = "#00B7EB", "SVM_RBF_LASSO" = "#E9C46A")

auc_data <- roc_data %>%
  group_by(BioProject, Model, ModelDisplay, Phenotype_Counts) %>%
  summarise(ROC_AUC = mean(ROC_AUC, na.rm = TRUE), .groups = "drop") %>%
  group_by(BioProject) %>%
  arrange(desc(ROC_AUC)) %>%
  mutate(rank = row_number(), ROC_AUC_Percent = ROC_AUC * 100,
         auc_label = sprintf("%s (%.1f%%)", ModelDisplay, ROC_AUC_Percent),
         label_x = 0.39, label_y = 0.38 - (rank - 1) * 0.085) %>%
  ungroup()

box_data <- auc_data %>% group_by(BioProject, Phenotype_Counts) %>% slice(1) %>% ungroup()

fig5B <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey55", linewidth = 0.55) +
  geom_line(linewidth = 1.4, alpha = 0.92) +
  facet_wrap(~ Phenotype_Counts, scales = "free", ncol = 4) +
  geom_rect(data = box_data, aes(xmin = 0.36, xmax = 1.00, ymin = 0.00, ymax = 0.44),
            inherit.aes = FALSE, fill = "white", color = "#bbbbbb", linewidth = 0.5) +
  geom_text(data = auc_data, aes(x = label_x, y = label_y, label = auc_label, color = Model),
            size = 3.0, fontface = "bold", hjust = 0, show.legend = FALSE) +
  scale_color_manual(values = model_colors,
                     labels = c("SVM_RBF" = "SVM-RBF (All Features)", "SVM_RBF_RFECV" = "SVM-RBF (RFECV)",
                                "SVM_RBF_PFI" = "SVM-RBF (PFI)", "SVM_RBF_LASSO" = "SVM-RBF (LASSO)")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01), labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01), labels = function(x) sprintf("%.1f", x)) +
  labs(title = "ROC Curve Comparison of SVM-RBF Models on Validation Datasets by BioProject",
       x = "False Positive Rate", y = "True Positive Rate", color = "Model") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "#1a1a1a", margin = margin(b = 8), lineheight = 1.3),
    axis.title.x = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 8)),
    axis.text = element_text(size = 10, face = "bold", color = "#000000"),
    axis.line = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.12, "cm"),
    strip.text = element_text(size = 10, face = "bold", color = "#000000", margin = margin(t = 5, b = 5), lineheight = 1.2),
    strip.background = element_rect(fill = "white", color = "#000000", linewidth = 0.7),
    panel.spacing = unit(0.9, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 8)),
    legend.text = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    legend.key.width = unit(1.4, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 6, b = 4),
    plot.margin = margin(12, 12, 10, 12)
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 2.5)))

print(fig5B)
ggsave(file.path(dir_main, "Fig5B_roc_bioproject_svm_rbf_validation.png"), fig5B, width = 16, height = 6, dpi = 600, bg = "white")
ggsave(file.path(dir_main, "Fig5B_roc_bioproject_svm_rbf_validation.pdf"), fig5B, width = 16, height = 6, device = "pdf", bg = "white")



# ==============================================================================
# SECTION 2 : SUPPLEMENTARY FIGURES
# ==============================================================================

# ------------------------------------------------------------------------------
# Supplementary Figure 1A - Classified vs Unclassified Reads per BioProject
# ------------------------------------------------------------------------------
df_reads <- read_delim("read_count_training_datasets_with_metadata.csv", delim = "\t")

df_reads <- df_reads %>%
  mutate(Total_Reads = Classified + Unclassified, Pct_Classified = (Classified / Total_Reads) * 100)

bio_avg <- df_reads %>%
  group_by(BioProject) %>%
  summarise(
    avg_classified = mean(Classified, na.rm = TRUE),
    avg_unclassified = mean(Unclassified, na.rm = TRUE),
    avg_total = mean(Classified + Unclassified, na.rm = TRUE),
    total_classified = sum(Classified, na.rm = TRUE),
    total_reads = sum(Classified + Unclassified, na.rm = TRUE),
    count = n(), .groups = "drop"
  ) %>%
  mutate(BioProject_label = paste0(BioProject, " (n=", count, ")"),
         pct_classified = (total_classified / total_reads) * 100) %>%
  arrange(avg_total)

bio_avg_long <- bio_avg %>%
  select(BioProject_label, avg_classified, avg_unclassified, avg_total, pct_classified) %>%
  pivot_longer(cols = c(avg_classified, avg_unclassified), names_to = "Category", values_to = "avg_reads") %>%
  mutate(Category = factor(Category, levels = c("avg_unclassified", "avg_classified"),
                           labels = c("Unclassified reads", "Classified reads")))

bio_avg_long$BioProject_label <- factor(bio_avg_long$BioProject_label, levels = unique(bio_avg$BioProject_label))

figS1A <- ggplot(bio_avg_long, aes(x = avg_reads / 1e6, y = BioProject_label, fill = Category)) +
  geom_col(data = subset(bio_avg_long, Category == "Classified reads"), position = "identity", width = 0.7, alpha = 0.85) +
  geom_col(data = subset(bio_avg_long, Category == "Unclassified reads"), position = "identity", width = 0.7, alpha = 0.85) +
  geom_text(data = bio_avg, aes(x = avg_total / 1e6, y = BioProject_label, label = paste0(round(pct_classified, 1), "%")),
            inherit.aes = FALSE, hjust = -0.1, color = "black", size = 3, fontface = "bold") +
  scale_fill_manual(name = "Read Type", values = c("Unclassified reads" = "#E57373", "Classified reads" = "#66BB6A")) +
  scale_x_continuous(labels = label_number(suffix = "M", scale = 1), expand = expansion(mult = c(0, 0.1)), breaks = pretty_breaks(n = 6)) +
  labs(
    title = "Average Classified and Unclassified Reads per BioProject",
    subtitle = "Each bar shows the mean number of reads per sample within each BioProject;\npercentages indicate classification rate",
    x = "Average number of reads per sample (millions)",
    y = "BioProject ID (n = number of samples)",
    caption = "Classified reads are taxonomically assigned; Unclassified reads lack taxonomic assignment"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5), color = "black"),
    plot.subtitle = element_text(size = 10, hjust = 0.5, margin = margin(b = 15), color = "grey30", lineheight = 1.2),
    plot.caption = element_text(size = 9, hjust = 0, color = "grey40", margin = margin(t = 10), face = "italic"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.4),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black", size = 9),
    axis.text.y = element_text(size = 8.5, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 12), color = "black"),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 12), color = "black"),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 11, face = "bold", color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.5),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.7, "cm"),
    legend.margin = margin(8, 10, 8, 10),
    legend.box.margin = margin(0, 0, 10, 0),
    plot.margin = margin(15, 30, 15, 15)
  )

print(figS1A)
ggsave(file.path(dir_supp, "FigS1A_bioproject_reads_publication.png"), figS1A, width = 13, height = 10, dpi = 600, bg = "white")


# ------------------------------------------------------------------------------
# Supplementary Figure 1B - Species Diversity vs Classified Reads Correlation
# ------------------------------------------------------------------------------
pearson_test  <- cor.test(df_reads$Species, df_reads$Classified, method = "pearson")
spearman_test <- cor.test(df_reads$Species, df_reads$Classified, method = "spearman")
lm_model      <- lm(Classified ~ Species, data = df_reads)
lm_summary    <- summary(lm_model)

figS1B <- ggplot(df_reads, aes(x = Species, y = Classified)) +
  geom_point(color = "#7570b3", size = 2.2, alpha = 0.5, shape = 16) +
  geom_smooth(method = "lm", aes(color = "Linear (Pearson)"), se = TRUE, fill = "#1b9e77", alpha = 0.15, linewidth = 1.1) +
  geom_smooth(method = "loess", aes(color = "Non-linear (Spearman)"), se = FALSE, linewidth = 1.1, linetype = "dashed") +
  scale_color_manual(
    name = "Regression Type",
    values = c("Linear (Pearson)" = "#1b9e77", "Non-linear (Spearman)" = "#d95f02"),
    labels = c(sprintf("Linear (Pearson r = %.3f***)", pearson_test$estimate),
               sprintf("Non-linear (Spearman \u03c1 = %.3f***)", spearman_test$estimate))
  ) +
  scale_y_continuous(labels = label_number(suffix = "M", scale = 1e-6, accuracy = 0.1)) +
  scale_x_continuous(labels = label_comma()) +
  labs(x = "Number of Species per Sample", y = "Classified Reads (Millions)",
       title = "Correlation Between Species Diversity and Classified Reads",
       subtitle = sprintf("R\u00b2 = %.3f, n = %d samples", lm_summary$r.squared, nrow(df_reads))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey20", margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.margin = margin(5, 10, 5, 10),
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(15, 15, 15, 15)
  )

print(figS1B)
ggsave(file.path(dir_supp, "FigS1B_species_vs_classified_reads.png"), figS1B, width = 10, height = 8, dpi = 600, bg = "white")

# Supporting per-BioProject read/species statistics (exported for Supplementary Table)
stats_classified <- df_reads %>%
  group_by(BioProject) %>%
  summarise(
    Mean_Classified = mean(Classified, na.rm = TRUE), SD_Classified = sd(Classified, na.rm = TRUE),
    Median_Classified = median(Classified, na.rm = TRUE), Min_Classified = min(Classified, na.rm = TRUE),
    Max_Classified = max(Classified, na.rm = TRUE), Mean_Species = mean(Species, na.rm = TRUE),
    SD_Species = sd(Species, na.rm = TRUE), Median_Species = median(Species, na.rm = TRUE),
    Sample_Count = n(), .groups = "drop"
  ) %>%
  arrange(desc(Sample_Count))

write.table(stats_classified, file.path(dir_supp, "SupplementaryTable_BioProject_Read_Species_Stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ------------------------------------------------------------------------------
# Supplementary Figure 6A - Species Retention Across Prevalence Filtering
# ------------------------------------------------------------------------------
df_prev <- read_tsv("Data_Metadata/otu_prevalence_filtering_data.tsv")
df_prev <- df_prev %>% filter(Prevalence_Threshold_Percent <= 25)

figS6A <- ggplot(df_prev, aes(x = Prevalence_Threshold_Percent, y = OTUs_Retained)) +
  geom_line(color = "#2C5F8D", linewidth = 1.2) +
  geom_point(size = 5, color = "#2C5F8D", fill = "white", shape = 21, stroke = 2) +
  geom_text(aes(label = scales::comma(OTUs_Retained)), vjust = -1.2, color = "#333333", fontface = "bold", size = 3.5) +
  labs(x = "Prevalence Threshold (% of Samples)", y = "Number of Species Retained",
       title = "Species Retention Across Prevalence Filtering Thresholds",
       subtitle = "Pooled analysis of 7,452 samples from 32 independent studies",
       caption = "Each point represents the number of species present in at least X% of samples") +
  scale_x_continuous(breaks = df_prev$Prevalence_Threshold_Percent, labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.05, 0.15))) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 15, color = "#1a1a1a", margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0, size = 11, color = "#4a4a4a", margin = margin(b = 15)),
    plot.caption = element_text(hjust = 0, size = 9, color = "#666666", face = "italic", margin = margin(t = 10)),
    axis.title.x = element_text(face = "bold", size = 13, color = "#000000", margin = margin(t = 12)),
    axis.title.y = element_text(face = "bold", size = 13, color = "#000000", margin = margin(r = 12)),
    axis.text.x = element_text(face = "bold", size = 11, color = "#000000", margin = margin(t = 5)),
    axis.text.y = element_text(face = "bold", size = 11, color = "#000000", margin = margin(r = 5)),
    panel.grid.major.y = element_line(color = "#e0e0e0", linewidth = 0.3),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    axis.line = element_line(color = "#1a1a1a", linewidth = 0.5),
    axis.ticks = element_line(color = "#1a1a1a", linewidth = 0.4),
    axis.ticks.length = unit(0.2, "cm"), panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(), plot.margin = margin(20, 20, 15, 15)
  )

print(figS6A)
ggsave(file.path(dir_supp, "FigS6A_prevalence_filtering_publication.png"), figS6A, width = 10, height = 6, dpi = 600, bg = "white")

# Species retention summary table for all prevalence thresholds
species_retention <- df_prev %>%
  select(Prevalence_Threshold_Percent, OTUs_Retained) %>%
  arrange(Prevalence_Threshold_Percent)
print(species_retention)


# ------------------------------------------------------------------------------
# Supplementary Figure 6B - F1 Score Across Prevalence Thresholds (0-25%)
# ------------------------------------------------------------------------------
metrics_file <- read.csv("perfomance_on_test_dataset_combined_for_all_models_prev_0_to_25.tsv", sep = "\t")

metrics_file <- metrics_file %>%
  filter(Prevalence <= 25, Feature_Selection == "All_Features", Model %in% c("SVM_RBF", "SVM_LIN", "RF", "LogReg")) %>%
  mutate(BaseModel = case_when(
    Model == "SVM_RBF" ~ "SVM-RBF", Model == "SVM_LIN" ~ "SVM-Linear",
    Model == "RF" ~ "RF", Model == "LogReg" ~ "LogReg-ElasticNet"
  ))

metrics_long <- metrics_file %>%
  pivot_longer(cols = "F1.Score", names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = "F1 Score", Value_Percent = Value * 100)

custom_colors <- c("SVM-RBF" = "#E63946", "SVM-Linear" = "#1D3557", "RF" = "#2A9D8F", "LogReg-ElasticNet" = "#F4A261")

max_points <- metrics_long %>% slice(which.max(Value_Percent))

figS6B <- ggplot(metrics_long, aes(x = Prevalence, y = Value_Percent, color = BaseModel, group = BaseModel)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3.5) +
  geom_point(data = max_points, aes(x = Prevalence, y = Value_Percent), color = "#000000", size = 6, shape = 8, stroke = 1.5) +
  geom_text(data = max_points, aes(x = Prevalence, y = Value_Percent, label = sprintf("%s (%.1f%%)", BaseModel, Value_Percent)),
            vjust = -1.5, color = "#000000", size = 4, fontface = "bold") +
  labs(title = "F1 Score Performance Across Prevalence Thresholds (0-25%)",
       subtitle = "Test Dataset: n = 1,491 (798 non-healthy, 693 healthy)",
       x = "Prevalence Threshold (%)", y = "F1 Score (%)", color = "Model") +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5), labels = paste0(seq(0, 25, by = 5), "%")) +
  scale_y_continuous(
    limits = c(floor(min(metrics_long$Value_Percent)) - 2, ceiling(max(metrics_long$Value_Percent)) + 3),
    breaks = seq(0, 100, by = 5), labels = function(x) paste0(x, "%")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0, color = "#1a1a1a", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 12, hjust = 0, color = "#1a1a1a", margin = margin(b = 12)),
    axis.title.x = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    axis.text = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line = element_line(color = "#000000", linewidth = 0.6),
    axis.ticks = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    panel.grid.major.y = element_line(color = "#e0e0e0", linewidth = 0.3),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.6),
    legend.position = "top", legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    legend.text = element_text(size = 12, face = "bold", margin = margin(r = 15)),
    legend.key.width = unit(1.8, "cm"), legend.key.height = unit(0.5, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 10, b = 0),
    plot.margin = margin(15, 15, 12, 15),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(figS6B)
ggsave(file.path(dir_supp, "FigS6B_F1_score_prev_0_25_all_models.png"), figS6B, width = 10, height = 6, dpi = 600, bg = "white")
ggsave(file.path(dir_supp, "FigS6B_F1_score_prev_0_25_all_models.pdf"), figS6B, width = 10, height = 6, device = "pdf", bg = "white")


# ------------------------------------------------------------------------------
# Supplementary Figure 6C - ROC Curve: SVM-RBF Top-50 Features (Test Dataset)
# ------------------------------------------------------------------------------
roc_data <- read.csv("roc_curve_testset_top50_svm_rbf.csv")

auc_value <- trapz(roc_data$FPR, roc_data$TPR)
auc_percent <- auc_value * 100

box_df <- data.frame(xmin = 0.36, xmax = 1.00, ymin = 0.00, ymax = 0.15)
auc_df <- data.frame(x = 0.39, y = 0.075, auc_label = sprintf("SVM-RBF (Top-50 Features)  %.1f%%", auc_percent))

figS6C <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.55) +
  geom_line(color = "#E63946", linewidth = 1.6, alpha = 0.92) +
  geom_rect(data = box_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "white", color = "#bbbbbb", linewidth = 0.4) +
  geom_text(data = auc_df, aes(x = x, y = y, label = auc_label), inherit.aes = FALSE,
            color = "#E63946", size = 4.2, fontface = "bold", hjust = 0) +
  labs(title = "ROC Curve of SVM-RBF Trained on Top-50 Features on Test Dataset (n = 1,491)",
       x = "False Positive Rate", y = "True Positive Rate") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0.01, 0.01), labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0.01, 0.01), labels = function(x) sprintf("%.1f", x)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "#1a1a1a", margin = margin(b = 4), lineheight = 1.3),
    axis.title.x = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 8)),
    axis.text = element_text(size = 11, face = "bold", color = "#000000"),
    axis.line = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks = element_line(color = "#000000", linewidth = 0.4),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(14, 14, 12, 14)
  )

print(figS6C)
ggsave(file.path(dir_supp, "FigS6C_roc_curve_svm_rbf_top50_features.png"), figS6C, width = 8, height = 7, dpi = 600, bg = "white")
ggsave(file.path(dir_supp, "FigS6C_roc_curve_svm_rbf_top50_features.pdf"), figS6C, width = 8, height = 7, device = "pdf", bg = "white")


# ------------------------------------------------------------------------------
# Supplementary Figure 6D - SVM-RBF Performance Heatmaps (Validation, per BioProject)
# ------------------------------------------------------------------------------
metrics_df <- read.csv("model_evaluation_on_each_bioproject_validation_datasets_20pct_2026.csv")

metrics_df <- subset(metrics_df, Samples >= 10)
metrics_df <- metrics_df %>%
  mutate(BioProject_Phenotype_Counts = str_replace(BioProject_Phenotype_Counts, ": ", "\n")) %>%
  mutate(AUC.Score = ifelse(is.na(AUC.Score) | !is.finite(AUC.Score), 0, AUC.Score)) %>%
  filter(grepl("^SVM_RBF", Model), !is.na(Balanced.Accuracy))

feature_counts <- c("All Features" = 4022, "RFECV" = 3472, "PFI" = 3127, "LASSO" = 828)

metrics_df <- metrics_df %>%
  mutate(
    Model_Label = case_when(
      Model == "SVM_RBF" ~ "All Features", grepl("_RFE", Model) ~ "RFECV",
      grepl("_PFI", Model) ~ "PFI", grepl("_LASSO", Model) ~ "LASSO", TRUE ~ "Other"
    ),
    Model_Label_Annotated = paste0(Model_Label, "\n(n=", feature_counts[Model_Label], ")"),
    Model_Label_Annotated = factor(Model_Label_Annotated,
                                   levels = paste0(c("All Features", "RFECV", "PFI", "LASSO"), "\n(n=", feature_counts[c("All Features", "RFECV", "PFI", "LASSO")], ")"))
  )

heatmap_data_long <- metrics_df %>%
  pivot_longer(cols = c(Accuracy, Balanced.Accuracy, Precision, Recall, F1.Score, AUC.Score),
               names_to = "Metric", values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(Value_Percent = Value * 100)

# Helper: derive sensible figure dimensions from data
get_dims <- function(data, is_faceted = FALSE) {
  n_rows <- length(unique(data$BioProject_Phenotype_Counts))
  n_cols <- length(unique(data$Model_Label_Annotated))
  n_metrics <- if (is_faceted) length(unique(data$Metric)) else 1
  list(w = max(7, n_cols * 1.8 + 2), h = max(5, n_rows * 0.55 * n_metrics + 2.5))
}

create_heatmap <- function(data, metric_filter = NULL, title, filename = NULL, filter_zero = TRUE) {
  if (!is.null(metric_filter)) data <- data %>% filter(Metric == metric_filter)
  if (filter_zero && !is.null(metric_filter) && metric_filter != "Accuracy") data <- data %>% filter(Value_Percent > 0)
  
  data <- data %>% mutate(text_color = ifelse(Value_Percent >= 65, "white", "#1a1a1a"))
  is_faceted <- is.null(metric_filter)
  dims <- get_dims(data, is_faceted)
  
  p <- ggplot(data, aes(x = Model_Label_Annotated, y = BioProject_Phenotype_Counts, fill = Value_Percent)) +
    geom_tile(color = "#000000", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.1f%%", Value_Percent), color = text_color), size = 3.8, fontface = "bold") +
    scale_fill_gradientn(colours = brewer.pal(9, "RdYlBu"), limits = c(0, 100), name = "Score (%)",
                         breaks = c(0, 25, 50, 75, 100), labels = c("0%", "25%", "50%", "75%", "100%")) +
    scale_color_identity() +
    labs(title = title, x = NULL, y = NULL)
  
  if (is_faceted) p <- p + facet_grid(Metric ~ ., scales = "free_y", space = "free_y")
  
  p <- p +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "#1a1a1a", margin = margin(b = 6), lineheight = 1.25),
      axis.text.x = element_text(size = 11, face = "bold", color = "#000000", angle = 0, hjust = 0.5, lineheight = 1.2),
      axis.text.y = element_text(size = 11, face = "bold", color = "#000000", lineheight = 1.1),
      axis.line = element_line(color = "#000000", linewidth = 0.5),
      axis.ticks = element_line(color = "#000000", linewidth = 0.4),
      legend.title = element_text(size = 10, face = "bold", color = "#000000", margin = margin(b = 6)),
      legend.text = element_text(size = 10, face = "bold", color = "#000000"),
      legend.key.height = unit(1.8, "cm"), legend.key.width = unit(0.55, "cm"),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.3),
      legend.margin = margin(l = 10),
      strip.text = element_text(size = 10, face = "bold", color = "#000000", margin = margin(t = 4, b = 4)),
      strip.background = element_rect(fill = "#f0f0f0", color = "#000000", linewidth = 0.4),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.spacing.y = unit(0.6, "lines"),
      plot.margin = margin(14, 15, 12, 14)
    )
  
  print(p)
  
  if (!is.null(filename)) {
    ggsave(file.path(dir_supp, paste0(filename, ".png")), p, width = dims$w, height = dims$h, dpi = 600, bg = "white")
    ggsave(file.path(dir_supp, paste0(filename, ".pdf")), p, width = dims$w, height = dims$h, device = "pdf", bg = "white")
  }
}

create_heatmap(heatmap_data_long, title = "SVM-RBF Performance Across All Metrics on Validation Datasets",
               filename = "FigS6D_validation_heatmap_all_metrics", filter_zero = FALSE)
create_heatmap(heatmap_data_long, metric_filter = "Accuracy", title = "SVM-RBF Accuracy on Validation Datasets",
               filename = "FigS6D_validation_heatmap_accuracy", filter_zero = FALSE)
create_heatmap(heatmap_data_long, metric_filter = "Balanced.Accuracy", title = "SVM-RBF Balanced Accuracy on Validation Datasets",
               filename = "FigS6D_validation_heatmap_balanced_accuracy")
create_heatmap(heatmap_data_long, metric_filter = "F1.Score", title = "SVM-RBF F1 Score on Validation Datasets",
               filename = "FigS6D_validation_heatmap_f1_score")
create_heatmap(heatmap_data_long, metric_filter = "AUC.Score", title = "SVM-RBF AUC Score on Validation Datasets",
               filename = "FigS6D_validation_heatmap_auc_score")

write.table(metrics_df, file.path(dir_supp, "SupplementaryTable_validation_heatmap_svm_rbf_data.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# End of script
# ==============================================================================
