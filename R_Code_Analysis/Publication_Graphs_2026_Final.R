setwd("C:/Users/bablu/OneDrive - Università degli Studi di Milano/CODE/Bablu_ML_Result")
library(pheatmap)
library(dplyr)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(forcats) 
library(RColorBrewer)
library(forcats) 
library(ggplot2)
library(readr)
library(ggplot2)
library(readr)


# Species Prevalence Filtering Analysis - Publication Figure

df <- read_tsv("Data_Metadata/otu_prevalence_filtering_data.tsv")
df <- df %>% filter(Prevalence_Threshold_Percent <= 25)

p_prevalence <- ggplot(df, aes(x = Prevalence_Threshold_Percent, y = OTUs_Retained)) +
  geom_line(color = "#2C5F8D", linewidth = 1.2) +
  geom_point(size = 5, color = "#2C5F8D", fill = "white", shape = 21, stroke = 2) +
  geom_text(aes(label = scales::comma(OTUs_Retained)), vjust = -1.2, color = "#333333", fontface = "bold", size = 3.5) +
  labs(x = "Prevalence Threshold (% of Samples)", y = "Number of Species Retained", 
       title = "Species Retention Across Prevalence Filtering Thresholds", 
       subtitle = "Pooled-analysis of 7,452 samples from 32 independent studies", 
       caption = "Each point represents the number of species present in at least X% of samples") +
  scale_x_continuous(breaks = df$Prevalence_Threshold_Percent, labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.05, 0.15))) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(plot.title = element_text(hjust = 0, face = "bold", size = 15, color = "#1a1a1a", margin = margin(b = 5)), 
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
        panel.border = element_blank(), plot.margin = margin(20, 20, 15, 15))
p_prevalence

ggsave("2026_Code_Publications/prevalence_filtering_publication.png", plot = p_prevalence, width = 10, height = 6, dpi = 600, bg = "white")
# Print species retention for all prevalence thresholds
species_retention <- df %>%
  select(Prevalence_Threshold_Percent, OTUs_Retained) %>%
  arrange(Prevalence_Threshold_Percent)

print("Species Retention across Prevalence Thresholds:")
print(species_retention)



# F1 Score Performance Across Prevalence Thresholds - Publication Figure

library(ggplot2)
library(dplyr)
library(tidyr)

metrics_file <- read.csv("perfomance_on_test_dataset_combined_for_all_models_prev_0_to_25.tsv", sep = "\t")

metrics_file <- metrics_file %>%
  filter(Prevalence <= 25, Feature_Selection == "All_Features", Model %in% c("SVM_RBF", "SVM_LIN", "RF", "LogReg")) %>%
  mutate(BaseModel = case_when(
    Model == "SVM_RBF" ~ "SVM-RBF",
    Model == "SVM_LIN" ~ "SVM-Linear",
    Model == "RF" ~ "RF",
    Model == "LogReg" ~ "LogReg-ElasticNet"
  ))

metrics_long <- metrics_file %>%
  pivot_longer(cols = "F1.Score", names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = "F1 Score", Value_Percent = Value * 100)

custom_colors <- c("SVM-RBF" = "#E63946", "SVM-Linear" = "#1D3557", "RF" = "#2A9D8F", "LogReg-ElasticNet" = "#F4A261")

main_plot <- ggplot(metrics_long, aes(x = Prevalence, y = Value_Percent, color = BaseModel, group = BaseModel)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3.5) +
  labs(title = "F1 Score Performance Across Prevalence Thresholds (0-25%)",
       subtitle = "Test Dataset: n = 1,491 (798 non-healthy, 693 healthy)",
       x = "Prevalence Threshold (%)", y = "F1 Score (%)", color = "Model") +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5), labels = paste0(seq(0, 25, by = 5), "%")) +
  scale_y_continuous(
    limits = c(floor(min(metrics_long$Value_Percent)) - 2, ceiling(max(metrics_long$Value_Percent)) + 3),
    breaks = seq(0, 100, by = 5),
    labels = function(x) paste0(x, "%")
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
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#000000", fill = NA, linewidth = 0.6),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    legend.text = element_text(size = 12, face = "bold", margin = margin(r = 15)),
    legend.key.width = unit(1.8, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 10, b = 0),
    plot.margin = margin(15, 15, 12, 15),
    panel.background = element_rect(fill = "white", color = NA)
  )

max_points <- metrics_long %>% slice(which.max(Value_Percent))

main_plot <- main_plot +
  geom_point(data = max_points, aes(x = Prevalence, y = Value_Percent),
             color = "#000000", size = 6, shape = 8, stroke = 1.5) +
  geom_text(data = max_points,
            aes(x = Prevalence, y = Value_Percent,
                label = sprintf("%s (%.1f%%)", BaseModel, Value_Percent)),
            vjust = -1.5, color = "#000000", size = 4, fontface = "bold")

print(main_plot)

# Save as PNG (high-resolution raster)
ggsave("2026_Code_Publications/ML_Pictures/F1_score_prev_0_25_SVM_RBF.png",
       plot = main_plot, width = 10, height = 6, dpi = 600, bg = "white")

ggsave("2026_Code_Publications/ML_Pictures/F1_score_prev_0_25_SVM_RBF.pdf",
       plot = main_plot, width = 10, height = 6, device = "pdf", bg = "white")



# F1 Score: Feature Selection Methods Performance - Publication Figure

library(ggplot2)
library(dplyr)

metrics_file <- read.csv("perfomance_on_test_dataset_combined_for_all_models_prev_0_to_25.tsv", sep = "\t")

metrics_file <- metrics_file %>%
  filter(Prevalence == 20, grepl("^(SVM_RBF|SVM_LIN|RF|LogReg)", Model)) %>%
  mutate(
    BaseModel = case_when(
      grepl("^SVM_RBF", Model) ~ "SVM-RBF",
      grepl("^SVM_LIN", Model) ~ "SVM-Linear",
      grepl("^RF", Model) ~ "RF",
      grepl("^LogReg", Model) ~ "LogReg-ElasticNet"
    ),
    FeatureSelection = case_when(
      Model %in% c("SVM_RBF", "SVM_LIN", "RF", "LogReg") ~ paste0("All Features\n(n=", Features, ")"),
      grepl("_LASSO$", Model) ~ paste0("LASSO\n(n=", Features, ")"),
      grepl("_RFECV$", Model) ~ paste0("RFECV\n(n=", Features, ")"),
      grepl("_PFI$", Model) ~ paste0("PFI\n(n=", Features, ")"),
      TRUE ~ paste0("All Features\n(n=", Features, ")")
    ),
    FeatureSelection = factor(FeatureSelection, levels = unique(FeatureSelection[order(Features, decreasing = TRUE)])),
    F1_Percent = F1.Score * 100
  )

custom_colors <- c("SVM-RBF" = "#E63946", "SVM-Linear" = "#1D3557", "RF" = "#2A9D8F", "LogReg-ElasticNet" = "#F4A261")

main_plot <- ggplot(metrics_file, aes(x = FeatureSelection, y = F1_Percent, fill = BaseModel)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75,
           color = "#000000", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%.1f%%", F1_Percent), group = BaseModel),
    position = position_dodge(width = 0.8),
    vjust = -0.6,          # lifted slightly above bar top
    size = 3.8,            # bigger than before (was 3.2)
    color = "#000000",
    fontface = "bold"
  ) +
  labs(
    title = "Model Preformance (F1 Score) on Test Dataset: n = 1,491 (798 non-healthy, 693 healthy)",
    x = "Feature Selection Method (Number of Features)",
    y = "F1 Score (%)",
    fill = "Model"
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18)),
    breaks = seq(0, 100, by = 10),
    limits = c(0, 100),
    labels = function(x) paste0(x, "%")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(size = 16, face = "bold", hjust = 0, color = "#1a1a1a", margin = margin(b = 3)),
    axis.title.x      = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 10)),
    axis.title.y      = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    axis.text.x       = element_text(size = 12, face = "bold", color = "#000000", angle = 0, hjust = 0.5),
    axis.text.y       = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    legend.position    = "top",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    legend.text        = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 12)),
    legend.key.size    = unit(0.5, "cm"),
    legend.key.width   = unit(1.3, "cm"),
    legend.background  = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin  = margin(t = 10, b = 0),
    plot.margin        = margin(12, 12, 10, 12),
    panel.background   = element_rect(fill = "white", color = NA)
  )

print(main_plot)

ggsave("2026_Code_Publications/ML_Pictures/F1_score_feature_selection_comparison.png",
       plot = main_plot, width = 11, height = 6, dpi = 600, bg = "white")

ggsave("2026_Code_Publications/ML_Pictures/F1_score_feature_selection_comparison.pdf",
       plot = main_plot, width = 11, height = 6, device = "pdf", bg = "white")




# Performance evaluation on Validation datasets;

library(ggplot2)
library(dplyr)

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

# Color palette
custom_colors <- c(
  "All Features" = "#E63946",
  "PFI"          = "#00B7EB",
  "RFECV"        = "#2A9D8F",
  "LASSO"        = "#E9C46A"
)

main_plot <- ggplot(metrics_file, aes(x = reorder(FeatureLabel, -F1_Percent), y = F1_Percent, fill = FeatureSelection)) +
  geom_bar(stat = "identity", width = 0.65, color = "#000000", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.1f%%", F1_Percent)),
    vjust  = -0.7,       # lifted clearly above bar
    size   = 5,          # bigger (was 4)
    color  = "#000000",
    fontface = "bold"
  ) +
  labs(
    title = "SVM-RBF F1 Score on Validation Dataset: n = 642 (426 non-healthy, 216 healthy)",
    x    = "Feature Selection Method",
    y    = "F1 Score (%)",
    fill = "Models"
  ) +
  scale_fill_manual(
    values = custom_colors,
    labels = c(
      "All Features" = "SVM-RBF (All Features)",
      "RFECV"        = "SVM-RBF (RFECV)",
      "PFI"          = "SVM-RBF (PFI)",
      "LASSO"        = "SVM-RBF (LASSO)"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18)),  # extra headroom so labels don't clip
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(size = 16, face = "bold", hjust = 0, color = "#1a1a1a", margin = margin(b = 4)),
    axis.title.y      = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    axis.title.x      = element_text(size = 12, face = "bold", color = "#000000", margin = margin(t = 10)),
    axis.text.x       = element_text(size = 12, face = "bold", color = "#000000"),
    axis.text.y       = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.border      = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    legend.position   = "top",
    legend.title      = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 10)),
    legend.text       = element_text(size = 12, face = "bold", color = "#000000", margin = margin(r = 12)),
    legend.key.width  = unit(1.5, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    panel.background  = element_rect(fill = "white", color = NA)
  )

print(main_plot)

ggsave("2026_Code_Publications/ML_Pictures/F1_score_validation_svm_rbf_all_methods.png",
       plot = main_plot, width = 12, height = 6, dpi = 600, bg = "white")

ggsave("2026_Code_Publications/ML_Pictures/F1_score_validation_svm_rbf_all_methods.pdf",
       plot = main_plot, width = 12, height = 6, device = "pdf", bg = "white")





# Health status: Class wise performance on Test datasets

library(dplyr)
library(tidyr)
library(ggplot2)

combined_df <- read.csv("classification_report_for_all_feature_selection_20pct_models.csv", sep = "\t")
combined_df <- combined_df %>%
  mutate(Prevalence = as.integer(gsub("pct", "", Prevalence)))

feature_counts <- c("All_Features" = 4022, "LASSO" = 828, "PFI" = 3127, "RFE" = 3472)

model_families <- list(
  SVM_RBF = c("SVM_RBF", "SVM_RBF_PFI", "SVM_RBF_RFECV", "SVM_RBF_LASSO"),
  SVM_LIN = c("SVM_LIN", "SVM_LIN_PFI", "SVM_LIN_RFECV", "SVM_LIN_LASSO"),
  RF      = c("RF", "RF_LASSO", "RF_PFI", "RF_RFECV"),
  LogReg  = c("LogReg", "LogReg_PFI", "LogReg_RFECV", "LogReg_LASSO")
)

family_colors <- c(SVM_RBF = "#E63946", SVM_LIN = "#1D3557", RF = "#2A9D8F", LogReg = "#F4A261")

model_to_family <- unlist(lapply(names(model_families), function(fam) {
  setNames(rep(fam, length(model_families[[fam]])), model_families[[fam]])
}))

df_combined <- combined_df %>%
  filter(Prevalence == 20, Dataset == "test", Class_Label %in% c("0", "1")) %>%
  mutate(
    Class_Label_Desc = case_when(
      Class_Label == "0" ~ "Non-Healthy",
      Class_Label == "1" ~ "Healthy",
      TRUE ~ Class_Label
    ),
    Model_Display = gsub("_RFECV", "_RFE", Model),
    Model_Display = gsub("SVM_RBF", "SVM-RBF", Model_Display),
    Model_Display = gsub("SVM_LIN", "SVM-Linear", Model_Display),
    Model_Display = gsub("_", " ", Model_Display),
    Model_Family = model_to_family[Model],
    F1_Percent = F1_Score * 100,
    Feature_Selection_Label = case_when(
      Feature_Selection == "All_Features" ~ paste0("All Features (n=", feature_counts["All_Features"], ")"),
      Feature_Selection == "LASSO"        ~ paste0("LASSO (n=",        feature_counts["LASSO"],        ")"),
      Feature_Selection == "PFI"          ~ paste0("PFI (n=",          feature_counts["PFI"],          ")"),
      Feature_Selection == "RFE"          ~ paste0("RFECV (n=",        feature_counts["RFE"],          ")")
    ),
    Feature_Selection_Label = factor(
      Feature_Selection_Label,
      levels = c(
        paste0("All Features (n=", feature_counts["All_Features"], ")"),
        paste0("LASSO (n=",        feature_counts["LASSO"],        ")"),
        paste0("PFI (n=",          feature_counts["PFI"],          ")"),
        paste0("RFECV (n=",        feature_counts["RFE"],          ")")
      )
    )
  ) %>%
  group_by(Feature_Selection_Label) %>%
  mutate(Model_Display = factor(
    Model_Display,
    levels = unique(Model_Display[order(Model_Family, Model_Display)])
  )) %>%
  ungroup()

p <- ggplot(df_combined, aes(x = Model_Display, y = F1_Percent, fill = Class_Label_Desc)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.7,
    color = "#000000",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", F1_Percent)),
    position = position_dodge(width = 0.75),
    vjust    = -0.7,
    size     = 4.4,     
    color    = "#000000",
    fontface = "bold"
  ) +
  facet_wrap(~ Feature_Selection_Label, nrow = 2, ncol = 2, scales = "free_x",) +
  scale_fill_manual(
    values = c("Non-Healthy" = "#E69F00", "Healthy" = "#56B4E9"),
    labels = c("Non-Healthy" = "Non-Healthy (n=798)", "Healthy" = "Healthy (n=693)")
  ) +
  labs(
    title    = "Class-Wise Models Perfromance (F1 Score) on Test Datasets (n=1,491)",
    x    = "",
    y    = "F1 Score (%)",
    fill = "Class"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.20)),  # slightly more headroom for bigger labels
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0, color = "#1a1a1a",
                                     margin = margin(b = 4), lineheight = 1.25),
    axis.title.y      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 10)),
    axis.text.x       = element_text(size = 11, face = "bold", color = "#000000",
                                     angle = 0, hjust = 0.5, vjust = 1),
    axis.text.y       = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    strip.text        = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(t = 5, b = 5)),
    strip.background  = element_rect(fill = "white", color = "#000000", linewidth = 0.7),
    panel.spacing     = unit(1.1, "lines"),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.border      = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.position   = "top",
    legend.title      = element_text(size = 14, face = "bold", color = "#000000",
                                     margin = margin(r = 10)),
    legend.text       = element_text(size = 14, face = "bold", color = "#000000",
                                     margin = margin(r = 12)),
    legend.key.width  = unit(2, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 5, b = 10),
    plot.margin       = margin(14, 14, 12, 14)
  )

print(p)

ggsave("2026_Code_Publications/ML_Pictures/f1_score_test_dataset_classwise_all_feature_selections.png",
       plot = p, width = 14, height = 8, dpi = 600, bg = "white")

ggsave("2026_Code_Publications/ML_Pictures/f1_score_test_dataset_classwise_all_feature_selections.pdf",
       plot = p, width = 14, height = 8, device = "pdf", bg = "white")



# Validation datasets Health status

library(dplyr)
library(tidyr)
library(ggplot2)

combined_df <- read.csv("classification_reports_validation_dataset_20pct_2026.tsv", sep = "\t")
combined_df <- combined_df %>%
  mutate(Prevalence = as.integer(gsub("pct", "", Prevalence)))

feature_counts <- c("All_Features" = 4022, "LASSO" = 828, "PFI" = 3127, "RFE" = 3472)
feature_selections_to_plot <- names(feature_counts)

df_filtered <- combined_df %>%
  filter(
    Prevalence == 20,
    Dataset == "validation",
    Class_Label %in% c("0", "1"),
    Feature_Selection %in% feature_selections_to_plot,
    Model %in% c("SVM_RBF", "SVM_RBF_PFI", "SVM_RBF_LASSO", "SVM_RBF_RFECV")
  ) %>%
  mutate(
    Class_Label_Desc = recode(Class_Label, "0" = "Non-Healthy", "1" = "Healthy"),
    Feature_Selection_Clean = recode(Feature_Selection,
                                     "RFE"          = "RFECV",
                                     "All_Features" = "All Features",
                                     .default = Feature_Selection),
    Feature_Selection_Annotated = paste0(Feature_Selection_Clean, "\n(n=", feature_counts[Feature_Selection], ")"),
    Feature_Selection_Annotated = factor(
      Feature_Selection_Annotated,
      levels = paste0(
        c("All Features", "RFECV", "PFI", "LASSO"),
        "\n(n=", feature_counts[c("All_Features", "RFE", "PFI", "LASSO")], ")"
      )
    ),
    F1_Percent = F1_Score * 100
  )

p_combined <- ggplot(df_filtered, aes(x = Feature_Selection_Annotated, y = F1_Percent, fill = Class_Label_Desc)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.7,
    color = "#000000",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", F1_Percent)),
    position = position_dodge(width = 0.75),
    vjust    = -0.7,     # lifted above bar top
    size     = 4.4,      # bigger, consistent with other figures
    color    = "#000000",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("Non-Healthy" = "#E69F00", "Healthy" = "#56B4E9"),
    labels = c("Non-Healthy" = "Non-Healthy (n=472)", "Healthy" = "Healthy (n=165)")
  ) +
  labs(
    x        = "",
    y        = "F1 Score (%)",
    fill     = "Class",
    title    = "Class-Wise SVM-RBF Models Performance (F1 Score) on Validation Dataset (n=642)",
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.20)),  # headroom for bigger labels
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5, color = "#1a1a1a",
                                     margin = margin(b = 4), lineheight = 1.3),
    axis.title.y      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 10)),
    axis.title.x      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(t = 10)),
    axis.text.x       = element_text(size = 12, face = "bold", color = "#000000"),
    axis.text.y       = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.border      = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.position   = "top",
    legend.title      = element_text(size = 14, face = "bold", color = "#000000",
                                     margin = margin(r = 10)),
    legend.text       = element_text(size = 14, face = "bold", color = "#000000",
                                     margin = margin(r = 12)),
    legend.key.width  = unit(1.5, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 5, b = 10),
    plot.margin       = margin(14, 14, 12, 14)
  )

print(p_combined)

# Save as PNG
ggsave("2026_Code_Publications/ML_Pictures/validation_f1_classwise_svm_rbf.png",
       plot = p_combined, width = 10, height = 6, dpi = 600, bg = "white")

# Save as PDF (vector)
ggsave("2026_Code_Publications/ML_Pictures/validation_f1_classwise_svm_rbf.pdf",
       plot = p_combined, width = 10, height = 6, device = "pdf", bg = "white")






# ROC Curves Across Feature Selection Methods - Publication Figure
library(ggplot2)
library(dplyr)
library(ggrepel)
library(pracma)

roc_data <- read.csv("roc_data_test_dataset_20pct.csv", stringsAsFactors = FALSE)


family_colors <- c("SVM_RBF" = "#E63946", "SVM_LIN" = "#1D3557", "RF" = "#2A9D8F", "LogReg" = "#F4A261")

roc_data <- roc_data %>%
  mutate(
    BaseModel = case_when(
      grepl("^SVM_RBF", Model) ~ "SVM_RBF",
      grepl("^SVM_LIN", Model) ~ "SVM_LIN",
      grepl("^RF",      Model) ~ "RF",
      grepl("^LogReg",  Model) ~ "LogReg",
      TRUE ~ Model
    ),
    Feature_Selection = case_when(
      Model %in% c("SVM_RBF", "SVM_LIN", "RF", "LogReg") ~ "All Features",
      grepl("_LASSO$", Model) ~ "LASSO",
      grepl("_RFECV$", Model) ~ "RFECV",
      grepl("_PFI$",   Model) ~ "PFI",
      TRUE ~ Feature_Selection
    ),
    Feature_Selection = factor(Feature_Selection,
                               levels = c("All Features", "LASSO", "PFI", "RFECV")),
    Feature_Selection_Label = case_when(
      Feature_Selection == "All Features" ~ "All Features (n=4,022)",
      Feature_Selection == "LASSO"        ~ "LASSO (n=828)",
      Feature_Selection == "PFI"          ~ "PFI (n=3,127)",
      Feature_Selection == "RFECV"        ~ "RFECV (n=3,472)"
    ),
    Feature_Selection_Label = factor(
      Feature_Selection_Label,
      levels = c("All Features (n=4,022)", "LASSO (n=828)",
                 "PFI (n=3,127)",          "RFECV (n=3,472)")
    ),
    DisplayModel = case_when(
      grepl("^SVM_RBF", Model) ~ "SVM-RBF",
      grepl("^SVM_LIN", Model) ~ "SVM-Linear",
      grepl("^RF",      Model) ~ "RF",
      grepl("^LogReg",  Model) ~ "LogReg (ElasticNet)"
    )
  )

# ── AUC table ─────────────────────────────────────────────────────────────────
auc_data <- roc_data %>%
  group_by(Model, DisplayModel, Feature_Selection, Feature_Selection_Label, BaseModel) %>%
  arrange(FPR) %>%
  summarise(AUC = trapz(FPR, TPR), .groups = "drop") %>%
  group_by(Feature_Selection_Label) %>%
  arrange(desc(AUC)) %>%
  mutate(
    rank        = row_number(),
    AUC_percent = AUC * 100,
    auc_label   = sprintf("%s (%.1f%%)", DisplayModel, AUC_percent),
    label_x     = 0.39,
    label_y     = 0.38 - (rank - 1) * 0.085
  ) %>%
  ungroup()

# Box extent — one row per panel
box_data <- auc_data %>%
  group_by(Feature_Selection_Label) %>%
  slice(1) %>%
  ungroup()

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = BaseModel)) +
  
  # Reference diagonal
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey55", linewidth = 0.55) +
  
  # ROC lines
  geom_line(linewidth = 1.4, alpha = 0.92) +
  
  # Facets — single row, no extra space
  facet_grid(. ~ Feature_Selection_Label) +
  
  # White background box for AUC table — bottom-right corner
  geom_rect(
    data        = box_data,
    aes(xmin = 0.36, xmax = 1.00, ymin = 0.00, ymax = 0.41),
    inherit.aes = FALSE,
    fill = "white", color = "#bbbbbb", linewidth = 0.5
  ) +
  
  # Combined Model (AUC%) label — single left-aligned call, no gap
  geom_text(
    data     = auc_data,aes(x = label_x, y = label_y, label = auc_label, color = BaseModel),
    size     = 4, fontface = "bold", hjust = 0,show.legend = FALSE) +
  labs(title = "ROC Curve Comparison of All ML Models on Test Dataset (n = 1,491)",
    x     = "False Positive Rate",y     = "True Positive Rate",color = "Model") +
  scale_color_manual(
    values = family_colors,
    labels = c("SVM_RBF" = "SVM-RBF", "SVM_LIN" = "SVM-Linear",
               "RF"      = "RF",      "LogReg"   = "LogReg (ElasticNet)")
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01),
                     labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01),
                     labels = function(x) sprintf("%.1f", x)) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5,
                                     color = "#1a1a1a", margin = margin(b = 8),
                                     lineheight = 1.3),
    axis.title.x      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(t = 8)),
    axis.title.y      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 8)),
    axis.text         = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.12, "cm"),
    strip.text        = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(t = 5, b = 5)),
    strip.background  = element_rect(fill = "white", color = "#000000", linewidth = 0.7),
    panel.spacing     = unit(0.8, "lines"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.title      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 8)),
    legend.text       = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 10)),
    legend.key.width  = unit(1.4, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 6, b = 4),
    plot.margin       = margin(12, 12, 10, 12)   # tighter vertical margins
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 2.5)))

print(p)

ggsave("2026_Code_Publications/ML_Pictures/roc_curves_ontest_datasets.png",
       plot = p, width = 14, height = 6, dpi = 600, bg = "white")
ggsave("2026_Code_Publications/ML_Pictures/roc_curves_ontest_datasets.pdf",
       plot = p, width = 18, height = 6, device = "pdf", bg = "white")




# ROC Curve - Validation Dataset (SVM-RBF, All Feature Selections)
library(ggplot2)
library(dplyr)
library(pracma)

roc_data <- read.csv("roc_data_validation_dataset_20pct_2026.tsv", sep = "\t", stringsAsFactors = FALSE)

roc_data$Model <- gsub("LR[-_]ElasticNet", "LogReg", roc_data$Model)


model_colors <- c(
  "SVM_RBF"       = "#E63946",
  "SVM_RBF_RFE"   = "#2A9D8F",
  "SVM_RBF_PFI"   = "#00B7EB",
  "SVM_RBF_LASSO" = "#E9C46A"
)

roc_data <- roc_data %>%
  filter(grepl("^SVM_RBF", Model)) %>%
  mutate(
    Feature_Selection = case_when(
      grepl("_RFE$",   Model) ~ "RFECV",
      grepl("_PFI$",   Model) ~ "PFI",
      grepl("_LASSO$", Model) ~ "LASSO",
      Model == "SVM_RBF"      ~ "All Features",
      TRUE ~ "Other"
    ),
    Feature_Selection = factor(Feature_Selection,
                               levels = c("All Features", "RFECV", "PFI", "LASSO")),
    DisplayModel = case_when(
      grepl("_RFE$",   Model) ~ "SVM-RBF (RFECV)",
      grepl("_PFI$",   Model) ~ "SVM-RBF (PFI)",
      grepl("_LASSO$", Model) ~ "SVM-RBF (LASSO)",
      Model == "SVM_RBF"      ~ "SVM-RBF (All Features)",
      TRUE ~ Model
    )
  )

# ── AUC table ─────────────────────────────────────────────────────────────────
auc_data <- roc_data %>%
  group_by(Model, DisplayModel, Feature_Selection) %>%
  arrange(FPR) %>%
  summarise(AUC = trapz(FPR, TPR), .groups = "drop") %>%
  arrange(desc(AUC)) %>%
  mutate(
    rank        = row_number(),
    AUC_percent = AUC * 100,
    auc_label   = sprintf("%s (%.1f%%)", DisplayModel, AUC_percent),
    label_x     = 0.5,
    label_y     = 0.38 - (rank - 1) * 0.085
  )

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
  
  # Reference diagonal
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey55", linewidth = 0.55) +
  
  # ROC lines
  geom_line(linewidth = 1.4, alpha = 0.92) +
  
  # White background box for AUC table — bottom-right corner
  geom_rect(
    data        = auc_data[1, ],
    aes(xmin = 0.49, xmax = 1.00, ymin = 0.1, ymax = 0.41),
    inherit.aes = FALSE,
    fill = "white", color = "#bbbbbb", linewidth = 0.5
  ) +
  
  # Combined Model (AUC%) label
  geom_text(
    data     = auc_data,
    aes(x = label_x, y = label_y, label = auc_label, color = Model),
    size     = 4, fontface = "bold", hjust = 0,
    show.legend = FALSE
  ) +
  
  labs(
    title = "ROC Curve Comparison of SVM-RBF Models on Validation Dataset (n = 642)",
    x     = "False Positive Rate",
    y     = "True Positive Rate",
    color = "Model"
  ) +
  
  scale_color_manual(
    values = model_colors,
    labels = c(
      "SVM_RBF"       = "SVM-RBF (All Features)",
      "SVM_RBF_RFE"   = "SVM-RBF (RFECV)",
      "SVM_RBF_PFI"   = "SVM-RBF (PFI)",
      "SVM_RBF_LASSO" = "SVM-RBF (LASSO)"
    )
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01),
                     labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01),
                     labels = function(x) sprintf("%.1f", x)) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5,
                                     color = "#1a1a1a", margin = margin(b = 8),
                                     lineheight = 1.3),
    axis.title.x      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(t = 8)),
    axis.title.y      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 8)),
    axis.text         = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.title      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 8)),
    legend.text       = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 10)),
    legend.key.width  = unit(1.4, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 6, b = 4),
    plot.margin       = margin(12, 12, 10, 12)
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 2.5)))

print(p)

ggsave("2026_Code_Publications/ML_Pictures/roc_validation_svm_rbf.png",
       plot = p, width = 12, height = 8, dpi = 600, bg = "white")
ggsave("2026_Code_Publications/ML_Pictures/roc_validation_svm_rbf.pdf",
       plot = p, width = 12, height = 8, device = "pdf", bg = "white")




# SVM-RBF Performance Heatmaps on Validation Datasets

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)

metrics_df <- read.csv("model_evaluation_on_each_bioproject_validation_datasets_20pct_2026.csv")

metrics_df <- subset(metrics_df, Samples >= 10)
metrics_df <- metrics_df %>%
  mutate(BioProject_Phenotype_Counts = str_replace(BioProject_Phenotype_Counts, ": ", "\n")) %>%
  mutate(AUC.Score = ifelse(is.na(AUC.Score) | !is.finite(AUC.Score), 0, AUC.Score))

metrics_df <- metrics_df %>%
  filter(grepl("^SVM_RBF", Model)) %>%
  filter(!is.na(Balanced.Accuracy))

feature_counts <- c("All Features" = 4022, "RFECV" = 3472, "PFI" = 3127, "LASSO" = 828)

metrics_df <- metrics_df %>%
  mutate(
    Model_Group = "SVM-RBF",
    Model_Label = case_when(
      Model == "SVM_RBF"      ~ "All Features",
      grepl("_RFE",   Model)  ~ "RFECV",
      grepl("_PFI",   Model)  ~ "PFI",
      grepl("_LASSO", Model)  ~ "LASSO",
      TRUE ~ "Other"
    ),
    Model_Label_Annotated = paste0(Model_Label, "\n(n=", feature_counts[Model_Label], ")"),
    Model_Label_Annotated = factor(
      Model_Label_Annotated,
      levels = paste0(
        c("All Features", "RFECV", "PFI", "LASSO"),
        "\n(n=", feature_counts[c("All Features", "RFECV", "PFI", "LASSO")], ")"
      )
    )
  )

heatmap_data_long <- metrics_df %>%
  pivot_longer(
    cols      = c(Accuracy, Balanced.Accuracy, Precision, Recall, F1.Score, AUC.Score),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) %>%
  mutate(Value_Percent = Value * 100)

# ── Helper: derive sensible figure dimensions from data ───────────────────────
get_dims <- function(data, is_faceted = FALSE) {
  n_rows <- length(unique(data$BioProject_Phenotype_Counts))
  n_cols <- length(unique(data$Model_Label_Annotated))
  n_metrics <- if (is_faceted) length(unique(data$Metric)) else 1
  w <- max(7,  n_cols * 1.8 + 2)
  h <- max(5,  n_rows * 0.55 * n_metrics + 2.5)
  list(w = w, h = h)
}

# ── Core plot function ────────────────────────────────────────────────────────
create_heatmap <- function(data, metric_filter = NULL, title,
                           filename = NULL, filter_zero = TRUE) {
  
  if (!is.null(metric_filter)) {
    data <- data %>% filter(Metric == metric_filter)
  }
  
  if (filter_zero && !is.null(metric_filter) && metric_filter != "Accuracy") {
    data <- data %>% filter(Value_Percent > 0)
  }
  
  # Text colour: white on dark tiles, black on light tiles
  data <- data %>%
    mutate(text_color = ifelse(Value_Percent >= 65, "white", "#1a1a1a"))
  
  is_faceted <- is.null(metric_filter)
  dims <- get_dims(data, is_faceted)
  
  p <- ggplot(data, aes(x = Model_Label_Annotated,
                        y = BioProject_Phenotype_Counts,
                        fill = Value_Percent)) +
    geom_tile(color = "#000000", linewidth = 0.35) +
    geom_text(
      aes(label = sprintf("%.1f%%", Value_Percent), color = text_color),
      size     = 3.8,        # bigger than original 3
      fontface = "bold"
    ) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "RdYlBu"),
      limits  = c(0, 100),
      name    = "Score (%)",
      breaks  = c(0, 25, 50, 75, 100),
      labels  = c("0%", "25%", "50%", "75%", "100%")
    ) +
    scale_color_identity() +
    labs(
      title = title,
      x     = NULL,
      y     = NULL
    )
  
  if (is_faceted) {
    p <- p + facet_grid(Metric ~ ., scales = "free_y", space = "free_y")
  }
  
  p <- p +
    theme_minimal(base_size = 12) +
    theme(
      plot.title        = element_text(size = 14, face = "bold", hjust = 0.5,
                                       color = "#1a1a1a", margin = margin(b = 6),
                                       lineheight = 1.25),
      axis.title.x      = element_text(size = 12, face = "bold", color = "#000000",
                                       margin = margin(t = 10)),
      axis.title.y      = element_text(size = 12, face = "bold", color = "#000000",
                                       margin = margin(r = 10)),
      axis.text.x       = element_text(size = 11, face = "bold", color = "#000000",
                                       angle = 0, hjust = 0.5, lineheight = 1.2),
      axis.text.y       = element_text(size = 11, face = "bold", color = "#000000",
                                       lineheight = 1.1),
      axis.line         = element_line(color = "#000000", linewidth = 0.5),
      axis.ticks        = element_line(color = "#000000", linewidth = 0.4),
      legend.title      = element_text(size = 10, face = "bold", color = "#000000",
                                       margin = margin(b = 6)),
      legend.text       = element_text(size = 10, face = "bold", color = "#000000"),
      legend.key.height = unit(1.8, "cm"),
      legend.key.width  = unit(0.55, "cm"),
      legend.position   = "right",
      legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.3),
      legend.margin     = margin(l = 10),
      strip.text        = element_text(size = 10, face = "bold", color = "#000000",
                                       margin = margin(t = 4, b = 4)),
      strip.background  = element_rect(fill = "#f0f0f0", color = "#000000", linewidth = 0.4),
      panel.grid        = element_blank(),
      panel.border      = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.spacing.y   = unit(0.6, "lines"),
      plot.margin       = margin(14, 15, 12, 14)
    )
  
  print(p)
  
  if (!is.null(filename)) {
    out_dir <- "2026_Code_Publications/ML_Pictures"
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    # PNG
    ggsave(file.path(out_dir, paste0(filename, ".png")),
           plot = p, width = dims$w, height = dims$h, dpi = 600, bg = "white")
    
    # PDF (vector)
    ggsave(file.path(out_dir, paste0(filename, ".pdf")),
           plot = p, width = dims$w, height = dims$h, device = "pdf", bg = "white")
  }
}

# ── Generate all heatmaps ─────────────────────────────────────────────────────

create_heatmap(
  heatmap_data_long,
  title       = "SVM-RBF Performance Across All Metrics on Validation Datasets",
  filename    = "validation_heatmap_svm_rbf_all_metrics",
  filter_zero = FALSE
)

create_heatmap(
  heatmap_data_long,
  metric_filter = "Accuracy",
  title         = "SVM-RBF Accuracy on Validation Datasets",
  filename      = "validation_heatmap_svm_rbf_accuracy",
  filter_zero   = FALSE
)

create_heatmap(
  heatmap_data_long,
  metric_filter = "Balanced.Accuracy",
  title         = "SVM-RBF Balanced Accuracy on Validation Datasets",
  filename      = "validation_heatmap_svm_rbf_balanced_accuracy"
)

create_heatmap(
  heatmap_data_long,
  metric_filter = "F1.Score",
  title         = "SVM-RBF F1 Score on Validation Datasets",
  filename      = "validation_heatmap_svm_rbf_f1_score"
)

create_heatmap(
  heatmap_data_long,
  metric_filter = "AUC.Score",
  title         = "SVM-RBF AUC Score on Validation Datasets",
  filename      = "validation_heatmap_svm_rbf_auc_score"
)

# ── Supplementary data export ─────────────────────────────────────────────────
write.table(
  metrics_df,
  file      = "2026_Code_Publications/ML_Pictures/validation_supplementary_heatmap_svm_rbf_data.tsv",
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)







# ROC Curves by BioProject - Validation Datasets (SVM-RBF)
library(ggplot2)
library(dplyr)
library(stringr)
library(pracma)

roc_data <- read.delim("roc_data_bioproject_20pct_2026.tsv", sep = "\t", header = TRUE)


roc_data <- roc_data %>%
  filter(grepl("^SVM_RBF", Model)) %>%
  mutate(
    Phenotype_Counts = str_replace_all(Phenotype_Counts, ": ", ":"),
    FeatureSelection = case_when(
      grepl("_RFE",   Model) ~ "RFECV",
      grepl("_PFI",   Model) ~ "PFI",
      grepl("_LASSO", Model) ~ "LASSO",
      Model == "SVM_RBF"     ~ "All Features",
      TRUE ~ "Other"
    ),
    ModelDisplay = case_when(
      grepl("_RFE",   Model) ~ "SVM-RBF (RFECV)",
      grepl("_PFI",   Model) ~ "SVM-RBF (PFI)",
      grepl("_LASSO", Model) ~ "SVM-RBF (LASSO)",
      Model == "SVM_RBF"     ~ "SVM-RBF (All Features)",
      TRUE ~ Model
    )
  )

model_colors <- c(
  "SVM_RBF"       = "#E63946",
  "SVM_RBF_RFECV" = "#2A9D8F",
  "SVM_RBF_PFI"   = "#00B7EB",
  "SVM_RBF_LASSO" = "#E9C46A"
)

auc_data <- roc_data %>%
  group_by(BioProject, Model, ModelDisplay, Phenotype_Counts) %>%
  summarise(ROC_AUC = mean(ROC_AUC, na.rm = TRUE), .groups = "drop") %>%
  group_by(BioProject) %>%
  arrange(desc(ROC_AUC)) %>%
  mutate(
    rank            = row_number(),
    ROC_AUC_Percent = ROC_AUC * 100,
    # Combined label: Model (AUC%) — consistent with other figures
    auc_label       = sprintf("%s (%.1f%%)", ModelDisplay, ROC_AUC_Percent),
    label_x         = 0.39,
    label_y         = 0.38 - (rank - 1) * 0.085
  ) %>%
  ungroup()

box_data <- auc_data %>%
  group_by(BioProject, Phenotype_Counts) %>%
  slice(1) %>%
  ungroup()

p2 <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
  
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "grey55", linewidth = 0.55) +
  
  geom_line(linewidth = 1.4, alpha = 0.92) +
  
  facet_wrap(~ Phenotype_Counts, scales = "free", ncol = 4) +
  
  geom_rect(
    data        = box_data,
    aes(xmin = 0.36, xmax = 1.00, ymin = 0.00, ymax = 0.44),
    inherit.aes = FALSE,
    fill = "white", color = "#bbbbbb", linewidth = 0.5
  ) +
  
  geom_text(
    data     = auc_data,
    aes(x = label_x, y = label_y, label = auc_label, color = Model),
    size     = 3.0, fontface = "bold", hjust = 0,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = model_colors,
    labels = c(
      "SVM_RBF"       = "SVM-RBF (All Features)",
      "SVM_RBF_RFECV" = "SVM-RBF (RFECV)",
      "SVM_RBF_PFI"   = "SVM-RBF (PFI)",
      "SVM_RBF_LASSO" = "SVM-RBF (LASSO)"
    )
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01),
                     labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01),
                     labels = function(x) sprintf("%.1f", x)) +
  
  labs(
    title = "ROC Curve Comparison of SVM-RBF Models on Validation Datasets by BioProjects ",
    x     = "False Positive Rate",
    y     = "True Positive Rate",
    color = "Model"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5,
                                     color = "#1a1a1a", margin = margin(b = 8),
                                     lineheight = 1.3),
    axis.title.x      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(t = 8)),
    axis.title.y      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 8)),
    axis.text         = element_text(size = 10, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.12, "cm"),
    strip.text        = element_text(size = 10, face = "bold", color = "#000000",
                                     margin = margin(t = 5, b = 5), lineheight = 1.2),
    strip.background  = element_rect(fill = "white", color = "#000000", linewidth = 0.7),
    panel.spacing     = unit(0.9, "lines"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.title      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 8)),
    legend.text       = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 10)),
    legend.key.width  = unit(1.4, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.background = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin = margin(t = 6, b = 4),
    plot.margin       = margin(12, 12, 10, 12)
  ) +
  guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 2.5)))

print(p2)

ggsave("2026_Code_Publications/ML_Pictures/roc_bioproject_svm_rbf_validation.png",
       plot = p2, width = 16, height = 6, dpi = 600, bg = "white")
ggsave("2026_Code_Publications/ML_Pictures/roc_bioproject_svm_rbf_validation.pdf",
       plot = p2, width = 16, height = 6, device = "pdf", bg = "white")




# ROC-AUC Curve: SVM-RBF Top 50 Features
library(ggplot2)
library(pracma)

roc_data <- read.csv("roc_curve_testset_top50_svm_rbf.csv")

auc_value   <- trapz(roc_data$FPR, roc_data$TPR)
auc_percent <- auc_value * 100

# AUC box — bottom-right (consistent with other ROC figures)
box_df  <- data.frame(xmin = 0.36, xmax = 1.00, ymin = 0.00, ymax = 0.15)
auc_df  <- data.frame(
  x         = 0.39,
  y         = 0.075,
  auc_label = sprintf("SVM-RBF (Top-50 Features)  %.1f%%", auc_percent)
)

p <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  
  # Reference diagonal
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey55", linewidth = 0.55) +
  
  # ROC line
  geom_line(color = "#E63946", linewidth = 1.6, alpha = 0.92) +
  
  # White AUC box
  geom_rect(
    data        = box_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "white", color = "#bbbbbb", linewidth = 0.4
  ) +
  
  # AUC label inside box
  geom_text(
    data        = auc_df,
    aes(x = x, y = y, label = auc_label),
    inherit.aes = FALSE,
    color    = "#E63946",
    size     = 4.2,
    fontface = "bold",
    hjust    = 0
  ) +
  
  labs(
    title    = "ROC Curve of SVM-RBF Trained on Top-50 Features on Test Dataset (n = 1,491)",
    x        = "False Positive Rate",
    y        = "True Positive Rate"
  ) +
  
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0.01, 0.01),
                     labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0.01, 0.01),
                     labels = function(x) sprintf("%.1f", x)) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5,
                                     color = "#1a1a1a", margin = margin(b = 4),
                                     lineheight = 1.3),
    plot.subtitle     = element_text(size = 11, hjust = 0.5, color = "#4a4a4a",
                                     margin = margin(b = 12), lineheight = 1.3),
    axis.title.x      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(t = 8)),
    axis.title.y      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 8)),
    axis.text         = element_text(size = 11, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.4),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    plot.margin       = margin(14, 14, 12, 14)
  )

print(p)

# Save as PNG
ggsave("2026_Code_Publications/ML_Pictures/roc_curve_svm_rbf_top50_features.png",
       plot = p, width = 8, height = 7, dpi = 600, bg = "white")

# Save as PDF (vector)
ggsave("2026_Code_Publications/ML_Pictures/roc_curve_svm_rbf_top50_features.pdf",
       plot = p, width = 8, height = 7, device = "pdf", bg = "white")




# Top 50 Species with Permutation Feature Importance (PFI)
library(readr)
library(dplyr)
library(ggplot2)

pfi <- read_tsv("PERMUTATION_IMPORTANCE_BACTERIAL_bacterial_20pct.tsv")
pfi <- pfi %>% filter(Importance_Mean > 0)

pfi_top50 <- pfi %>%
  arrange(desc(Importance_Mean)) %>%
  slice_head(n = 50)

p <- ggplot(pfi_top50, aes(x = Importance_Mean,
                           y = reorder(Feature, Importance_Mean),
                           fill = Importance_Mean)) +
  geom_col(color = "#000000", width = 0.75, linewidth = 0.35) +
  geom_errorbar(
    aes(xmin = pmax(0, Importance_Mean - Importance_Std),
        xmax = Importance_Mean + Importance_Std),
    width = 0.4, linewidth = 0.5, color = "#000000"
  ) +
  geom_text(
    aes(label = sprintf("%.4f", Importance_Mean)),
    hjust    = -0.15,
    size     = 3.2,          # bigger than original 2.8
    color    = "#000000",
    fontface = "bold"
  ) +
  
  # 3-colour gradient: low = blue, mid = white, high = red
  scale_fill_gradientn(
    colors = c("#3182bd", "#FFFFFF", "#e34a33"),
    name   = "Importance Score"
  ) +
  
  labs(
    title = "Top 50 Bacterial Species by SVM-RBF-PFI",
    x     = "Mean Permutation Importance Score (± Std)",
    y     = ""
  ) +
  
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5,
                                     color = "#1a1a1a", margin = margin(b = 6),
                                     lineheight = 1.3),
    axis.title.x      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(t = 10)),
    axis.text.y       = element_text(size = 12, face = "bold.italic",
                                     color = "#000000"),   # italic for species names
    axis.text.x       = element_text(size = 11, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.4),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major.x = element_line(color = "white", linewidth = 0.25),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(14, 14, 12, 14),
    legend.position    = "top",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 11, face = "bold", color = "#000000",
                                      margin = margin(r = 8)),
    legend.text        = element_text(size = 10, face = "bold", color = "#000000"),
    legend.key.width   = unit(3.0, "cm"),
    legend.key.height  = unit(0.45, "cm"),
    legend.background  = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin  = margin(t = 5, b = 10)
  )

print(p)
print(pfi_top50, n = 50)

# Save as PNG
ggsave("2026_Code_Publications/ML_Pictures/permutation_importance_top50_svm_rbf_pfi.png",
       plot = p, width = 12, height = 14, dpi = 600, bg = "white")

# Save as PDF (vector)
ggsave("2026_Code_Publications/ML_Pictures/permutation_importance_top50_svm_rbf_pfi.pdf",
       plot = p, width = 12, height = 14, device = "pdf", bg = "white")





# Maaslin2

library(pheatmap)
library(dplyr)
library(tidyr)
library(grid)

sig_df <- read.delim("top_50_pfi_maaslin_significant.tsv", sep = "\t")

maaslin2_heatmap <- function(
    df,
    title = "Top 50 Microbial Features: Significant Associations\n(-log(qval) * sign(coef))",
    cell_value = 'qval',
    data_label = 'data',
    metadata_label = 'metadata',
    border_color = '#000000',
    color = colorRampPalette(c("#3182bd", "#FFFFFF", "#e34a33")),
    col_rotate = 90,
    first_n = 50
) {
  
  if (!is.na(first_n) & first_n > 0 & first_n < nrow(df)) {
    if (cell_value == 'coef') {
      df <- df[order(-abs(df[[cell_value]])), ]
    } else {
      df <- df[order(df[[cell_value]]), ]
    }
    df_sub <- df[1:first_n, ]
    for (first_n_index in seq(first_n, nrow(df))) {
      if (length(unique(df_sub$feature)) == first_n) {
        break
      }
      df_sub <- df[1:first_n_index, ]
    }
    df <- df[which(df$feature %in% df_sub$feature), ]
  }
  
  if (nrow(df) < 2) {
    message('There are no associations to plot!')
    return(NULL)
  }
  
  metadata <- df$metadata
  data <- df$feature
  dfvalue <- df$value
  value <- NA
  
  if (cell_value == "pval") {
    value <- -log(df$pval) * sign(df$coef)
    value <- pmax(-20, pmin(20, value))
  } else if (cell_value == "qval") {
    value <- -log(df$qval) * sign(df$coef)
    value <- pmax(-20, pmin(20, value))
  } else if (cell_value == "coef") {
    value <- df$coef
  }
  
  verbose_metadata <- c()
  metadata_multi_level <- c()
  for (i in unique(metadata)) {
    levels <- unique(df$value[df$metadata == i])
    if (length(levels) > 1) {
      metadata_multi_level <- c(metadata_multi_level, i)
      for (j in levels) {
        verbose_metadata <- c(verbose_metadata, j)
      }
    } else {
      if (i == "Health_status") {
        verbose_metadata <- c(verbose_metadata, "Healthy")
      } else {
        verbose_metadata <- c(verbose_metadata, i)
      }
    }
  }
  
  n <- length(unique(data))
  m <- length(unique(verbose_metadata))
  
  if (n < 2 || m < 2) {
    message('Not enough features or metadata to create a heatmap plot.')
    return(NULL)
  }
  
  a <- matrix(0, nrow = n, ncol = m)
  a <- as.data.frame(a)
  rownames(a) <- unique(data)
  colnames(a) <- unique(verbose_metadata)
  
  for (i in seq_len(nrow(df))) {
    current_metadata <- metadata[i]
    if (current_metadata %in% metadata_multi_level) {
      current_metadata <- dfvalue[i]
    } else if (current_metadata == "Health_status") {
      current_metadata <- "Healthy"
    }
    if (abs(a[as.character(data[i]), as.character(current_metadata)]) > abs(value[i])) {
      next
    }
    a[as.character(data[i]), as.character(current_metadata)] <- value[i]
  }
  
  max_value <- max(abs(a), na.rm = TRUE)
  if (max_value == 0) max_value <- 1
  breaks <- seq(-max_value, max_value, length.out = 101)
  color_gradient <- color(length(breaks))
  
  p <- pheatmap::pheatmap(
    a,
    cellwidth = 20,
    cellheight = 12,
    main = title,
    fontsize = 13,
    fontsize_row = 9,
    fontsize_col = 11,
    border = TRUE,
    border_color = border_color,
    show_rownames = TRUE,
    show_colnames = TRUE,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    legend = TRUE,
    legend_breaks = seq(-ceiling(max_value), ceiling(max_value), by = max(1, ceiling(max_value/5))),
    legend_labels = sprintf("%.1f", seq(-ceiling(max_value), ceiling(max_value), by = max(1, ceiling(max_value/5)))),
    color = color_gradient,
    breaks = breaks,
    treeheight_row = 50,
    treeheight_col = 50,
    display_numbers = matrix(ifelse(a > 0.0, "+", ifelse(a < 0.0, "-", "")), nrow(a)),
    number_color = "#000000",
    fontsize_number = 8,
    angle_col = "90",
    silent = FALSE
  )
  
  message("Heatmap value range: ", round(min(a, na.rm = TRUE), 2), " to ", round(max(a, na.rm = TRUE), 2))
  return(p)
}

heatmap_plot <- maaslin2_heatmap(
  df = sig_df,
  cell_value = "qval",
  border_color = "#000000",
  color = colorRampPalette(c("#3182bd", "#FFFFFF", "#e34a33")),
  first_n = 50
)

png("2026_Code_Publications/ML_Pictures/maaslin2_heatmap_top50_features.png", 
    width = 4000, height = 4500, res = 400, type = "cairo")

grid.newpage()
grid.draw(heatmap_plot$gtable)

dev.off()





#### Leave-One-Study-Out (LOSO) Cross-Validation - Accuracy
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

df <- read_tsv("Leaveoneout_performance_test.tsv")

df <- df %>%
  filter(!is.na(Accuracy_mean_CV)) %>%
  mutate(
    BioProject_Label  = str_replace(Summary, ": ", "\n"),
    Acc_test_pct      = Accuracy_test     * 100,
    Acc_mean_CV_pct   = Accuracy_mean_CV  * 100,
    Acc_std_CV_pct    = Accuracy_std_CV   * 100,
    CV_lower          = Acc_mean_CV_pct - Acc_std_CV_pct,
    CV_upper          = Acc_mean_CV_pct + Acc_std_CV_pct,
    LOSO_better       = Acc_test_pct > Acc_mean_CV_pct
  ) %>%
  mutate(BioProject_Label = factor(BioProject_Label, levels = BioProject_Label))

mean_test_pct <- mean(df$Acc_test_pct,    na.rm = TRUE)
sd_test_pct   <- sd(df$Acc_test_pct,      na.rm = TRUE)
mean_cv_pct   <- mean(df$Acc_mean_CV_pct, na.rm = TRUE)

p <- ggplot(df, aes(x = BioProject_Label)) +
  
  # ── CV band (±SD shading) ─────────────────────────────────────────────────
  geom_rect(
    aes(xmin = as.numeric(BioProject_Label) - 0.4,
        xmax = as.numeric(BioProject_Label) + 0.4,
        ymin = CV_lower, ymax = CV_upper),
    fill = "#1D3557", alpha = 0.12
  ) +
  
  # ── Overall mean reference lines ─────────────────────────────────────────
  geom_hline(yintercept = mean_test_pct,
             linetype = "dashed", color = "#E63946", linewidth = 0.9, alpha = 0.85) +
  geom_hline(yintercept = mean_cv_pct,
             linetype = "dashed", color = "#1D3557", linewidth = 0.9, alpha = 0.85) +
  
  # ── CV error bars ─────────────────────────────────────────────────────────
  geom_segment(
    aes(x    = as.numeric(BioProject_Label),
        xend = as.numeric(BioProject_Label),
        y    = CV_lower, yend = CV_upper),
    color = "#1D3557", linewidth = 1.1, alpha = 0.7
  ) +
  
  # ── CV mean point ─────────────────────────────────────────────────────────
  geom_point(
    aes(y = Acc_mean_CV_pct,
        fill  = "Within-Study CV (Mean ± SD)",
        shape = "Within-Study CV (Mean ± SD)"),
    size = 5, color = "#000000", stroke = 0.6, alpha = 0.95
  ) +
  
  # ── LOSO test point ───────────────────────────────────────────────────────
  geom_point(
    aes(y = Acc_test_pct,
        fill  = "LOSO Test Performance",
        shape = "LOSO Test Performance"),
    size = 5.5, color = "#000000", stroke = 0.6, alpha = 1
  ) +
  
  # ── Star label for studies where LOSO > CV ───────────────────────────────
  geom_text(
    data = subset(df, LOSO_better == TRUE),
    aes(y = Acc_test_pct,
        label = paste0("★\n", round(Acc_test_pct, 1), "%")),
    vjust = -0.5, hjust = 0.5,
    size = 3.5, fontface = "bold", color = "#E63946"
  ) +
  
  # ── Mean annotations ──────────────────────────────────────────────────────
  annotate("text",
           x = nrow(df) - 0.3, y = mean_test_pct + 3.2,
           label  = sprintf("Mean LOSO: %.1f%%", mean_test_pct),
           hjust  = 1, vjust = 0,
           size   = 4.2, fontface = "bold", color = "#E63946") +
  
  annotate("text",
           x = nrow(df) - 0.3, y = mean_cv_pct - 3.2,
           label  = sprintf("Mean CV: %.1f%%", mean_cv_pct),
           hjust  = 1, vjust = 1,
           size   = 4.2, fontface = "bold", color = "#1D3557") +
  
  # ── Scales ────────────────────────────────────────────────────────────────
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  
  scale_fill_manual(
    name   = "",
    values = c("Within-Study CV (Mean ± SD)" = "#1D3557",
               "LOSO Test Performance"        = "#E63946"),
    labels = c("Within-Study CV (Mean ± SD)" = "Within-Study CV (Mean ± SD)",
               "LOSO Test Performance"        = "LOSO Test — Generalisation to unseen study")
  ) +
  
  scale_shape_manual(
    name   = "",
    values = c("Within-Study CV (Mean ± SD)" = 21,
               "LOSO Test Performance"        = 24),
    labels = c("Within-Study CV (Mean ± SD)" = "Within-Study CV (Mean ± SD)",
               "LOSO Test Performance"        = "LOSO Test — Generalisation to unseen study")
  ) +
  
  labs(
    title    = "Leave-One-Study-Out (LOSO) Cross-Validation: Model Generalisation Across Studies",
    subtitle = "Each study held out once as the test set | Blue band: Within-study CV variability (±SD) | ★ LOSO > CV mean",
    x        = "",
    y        = "Accuracy (%)"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5,
                                     color = "#1a1a1a", margin = margin(b = 4),
                                     lineheight = 1.3),
    plot.subtitle     = element_text(size = 10, hjust = 0.5, color = "#1a1a1a",
                                     margin = margin(b = 14), lineheight = 1.3),
    axis.title.y      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 10)),
    axis.text.x       = element_text(size = 9, face = "bold", color = "#000000",
                                     angle = 45, hjust = 1, vjust = 1,
                                     lineheight = 1.1),
    axis.text.y       = element_text(size = 11, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    legend.position    = "top",
    legend.direction   = "horizontal",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 11, face = "bold", color = "#000000"),
    legend.key.size    = unit(0.65, "cm"),
    legend.background  = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin  = margin(t = 5, b = 10),
    plot.margin        = margin(14, 14, 12, 14)
  ) +
  
  guides(
    fill  = guide_legend(override.aes = list(size = 5.5, stroke = 0.5)),
    shape = guide_legend(override.aes = list(size = 5.5, stroke = 0.5))
  )

print(p)

cat(sprintf("\nTotal studies included: %d\n",          nrow(df)))
cat(sprintf("Mean LOSO Accuracy:           %.1f%%\n",  mean_test_pct))
cat(sprintf("Mean Within-Study CV Accuracy: %.1f%%\n", mean_cv_pct))

# Save as PNG
ggsave("2026_Code_Publications/ML_Pictures/loso_accuracy_all_studies.png",
       plot = p, width = 16, height = 8, dpi = 600, bg = "white")

# Save as PDF (vector)
ggsave("2026_Code_Publications/ML_Pictures/loso_accuracy_all_studies.pdf",
       plot = p, width = 16, height = 8, device = "pdf", bg = "white")




#### Leave-One-Study-Out (LOSO) Cross-Validation - F1 Score (Filtered_main_text)
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

df <- read_tsv("Leaveoneout_performance_test.tsv")

df <- df %>%
  filter(BioProject != "PRJEB18755") %>%
  filter(!is.na(F1_test)) %>%
  mutate(BioProject_Label = str_wrap(Summary, width = 25)) %>%
  mutate(BioProject_Label = factor(BioProject_Label, levels = BioProject_Label)) %>%
  mutate(
    F1_test_pct    = F1_test    * 100,
    F1_mean_CV_pct = F1_mean_CV * 100,
    F1_std_CV_pct  = F1_std_CV  * 100,
    CV_lower       = F1_mean_CV_pct - F1_std_CV_pct,
    CV_upper       = F1_mean_CV_pct + F1_std_CV_pct,
    LOSO_better    = F1_test_pct > F1_mean_CV_pct
  )

mean_test_pct <- mean(df$F1_test_pct,    na.rm = TRUE)
sd_test_pct   <- sd(df$F1_test_pct,      na.rm = TRUE)
mean_cv_pct   <- mean(df$F1_mean_CV_pct, na.rm = TRUE)

p <- ggplot(df, aes(x = BioProject_Label)) +
  
  # ── CV band (±SD shading) ─────────────────────────────────────────────────
  geom_rect(
    aes(xmin = as.numeric(BioProject_Label) - 0.4,
        xmax = as.numeric(BioProject_Label) + 0.4,
        ymin = CV_lower, ymax = CV_upper),
    fill = "#1D3557", alpha = 0.12
  ) +
  
  # ── Overall mean reference lines ─────────────────────────────────────────
  geom_hline(yintercept = mean_test_pct,
             linetype = "dashed", color = "#E63946", linewidth = 0.9, alpha = 0.85) +
  geom_hline(yintercept = mean_cv_pct,
             linetype = "dashed", color = "#1D3557", linewidth = 0.9, alpha = 0.85) +
  
  # ── CV error bars ─────────────────────────────────────────────────────────
  geom_segment(
    aes(x    = as.numeric(BioProject_Label),
        xend = as.numeric(BioProject_Label),
        y    = CV_lower, yend = CV_upper),
    color = "#1D3557", linewidth = 1.1, alpha = 0.7
  ) +
  
  # ── CV mean point ─────────────────────────────────────────────────────────
  geom_point(
    aes(y     = F1_mean_CV_pct,
        fill  = "Within-Study CV (Mean ± SD)",
        shape = "Within-Study CV (Mean ± SD)"),
    size = 5, color = "#000000", stroke = 0.6, alpha = 0.95
  ) +
  
  # ── LOSO test point ───────────────────────────────────────────────────────
  geom_point(
    aes(y     = F1_test_pct,
        fill  = "LOSO Test Performance",
        shape = "LOSO Test Performance"),
    size = 5.5, color = "#000000", stroke = 0.6, alpha = 1
  ) +
  
  # ── Star label for studies where LOSO > CV ───────────────────────────────
  geom_text(
    data = subset(df, LOSO_better == TRUE),
    aes(y = F1_test_pct,
        label = paste0("★\n", round(F1_test_pct, 1), "%")),
    vjust = -0.5, hjust = 0.5,
    size = 3.5, fontface = "bold", color = "#E63946"
  ) +
  
  # ── Mean annotations ──────────────────────────────────────────────────────
  annotate("text",
           x = nrow(df) - 0.3, y = mean_test_pct + 3.2,
           label  = sprintf("Mean LOSO F1: %.1f%%", mean_test_pct),
           hjust  = 1, vjust = 0,
           size   = 4.2, fontface = "bold", color = "#E63946") +
  
  annotate("text",
           x = nrow(df) - 0.3, y = mean_cv_pct - 3.2,
           label  = sprintf("Mean CV F1: %.1f%%", mean_cv_pct),
           hjust  = 1, vjust = 1,
           size   = 4.2, fontface = "bold", color = "#1D3557") +
  
  # ── Scales ────────────────────────────────────────────────────────────────
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  
  scale_fill_manual(
    name   = "",
    values = c("Within-Study CV (Mean ± SD)" = "#1D3557",
               "LOSO Test Performance"        = "#E63946"),
    labels = c("Within-Study CV (Mean ± SD)" = "Within-Study CV (Mean ± SD)",
               "LOSO Test Performance"        = "LOSO Test — Generalisation to unseen study")
  ) +
  
  scale_shape_manual(
    name   = "",
    values = c("Within-Study CV (Mean ± SD)" = 21,
               "LOSO Test Performance"        = 24),
    labels = c("Within-Study CV (Mean ± SD)" = "Within-Study CV (Mean ± SD)",
               "LOSO Test Performance"        = "LOSO Test — Generalisation to unseen study")
  ) +
  
  labs(
    title    = "Leave-One-Study-Out (LOSO) Cross-Validation: F1 Score Generalisation Across Studies",
    x        = "",
    y        = "F1 Score (%)"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5,
                                     color = "#1a1a1a", margin = margin(b = 4),
                                     lineheight = 1.3),

    axis.title.y      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 10)),
    axis.text.x       = element_text(size = 12, face = "bold", color = "#000000",
                                     angle = 45, hjust = 1, vjust = 1,
                                     lineheight = 1.1),
    axis.text.y       = element_text(size = 12, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    legend.position    = "top",
    legend.direction   = "horizontal",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 12, face = "bold", color = "#000000"),
    legend.key.size    = unit(0.65, "cm"),
    legend.background  = element_rect(fill = "white", color = "#000000", linewidth = 0.5),
    legend.box.margin  = margin(t = 5, b = 10),
    plot.margin        = margin(14, 14, 12, 14)
  ) +
  
  guides(
    fill  = guide_legend(override.aes = list(size = 5.5, stroke = 0.5)),
    shape = guide_legend(override.aes = list(size = 5.5, stroke = 0.5))
  )

print(p)

cat(sprintf("\nTotal studies included: %d\n",          nrow(df)))
cat(sprintf("Mean LOSO F1 Score:            %.1f%%\n", mean_test_pct))
cat(sprintf("Mean Within-Study CV F1 Score: %.1f%%\n", mean_cv_pct))

ggsave("2026_Code_Publications/ML_Pictures/loso_f1_score_filtered_studies.png",
       plot = p, width = 16, height = 8, dpi = 600, bg = "white")
ggsave("2026_Code_Publications/ML_Pictures/loso_f1_score_filtered_studies.pdf",
       plot = p, width = 16, height = 8, device = "pdf", bg = "white")






#### Leave-One-Study-Out ROC Curves - Publication Figure
library(ggplot2)
library(readr)
library(dplyr)
library(pracma)

roc_dir <- "LOB_ROC/"

bio_info <- list(
  "PRJEB27928"  = list(disease = "CRC", healthy = 119, disease_n = 139),
  "PRJEB53401"  = list(disease = "PD",  healthy =  70, disease_n =  48),
  "PRJEB6337"   = list(disease = "LV",  healthy = 142, disease_n = 163),
  "PRJNA375935" = list(disease = "AS",  healthy = 114, disease_n =  94),
  "PRJNA743718" = list(disease = "PD",  healthy =  18, disease_n =   8),
  "PRJNA893901" = list(disease = "CD n=31 / UC n=32",
                       cd_n = 31, uc_n = 32, healthy = 42, disease_n = 63)
)

roc_files      <- names(bio_info)
roc_file_names <- paste0("ROC_", roc_files, ".tsv")

# Distinct, colour-blind-friendly palette for 6 studies
colors <- c(
  "PRJEB27928"  = "#E63946",
  "PRJEB53401"  = "#1D3557",
  "PRJEB6337"   = "#2A9D8F",
  "PRJNA375935" = "#F4A261",
  "PRJNA743718" = "#9B59B6",
  "PRJNA893901" = "#457B9D"
)

roc_data_all <- data.frame()
auc_labels   <- data.frame()

for (i in seq_along(roc_file_names)) {
  f   <- roc_file_names[i]
  df  <- read_tsv(file.path(roc_dir, f), show_col_types = FALSE)
  bp  <- roc_files[i]
  info <- bio_info[[bp]]
  
  # Normalise to 0-1 if needed
  if (max(df$FPR, na.rm = TRUE) > 1) df$FPR <- df$FPR / 100
  if (max(df$TPR, na.rm = TRUE) > 1) df$TPR <- df$TPR / 100
  
  df <- df %>% arrange(FPR)
  auc_val <- trapz(df$FPR, df$TPR)
  
  df$BioProject <- bp
  df$Disease    <- info$disease
  roc_data_all  <- rbind(roc_data_all, df)
  
  auc_labels <- rbind(auc_labels, data.frame(
    BioProject = bp,
    Disease    = info$disease,
    AUC_pct    = auc_val * 100,
    # Clean label: BioProject · Disease (n_disease / n_healthy) · AUC
    Full_Label = sprintf("%s (%s=%d, Healthy=%d)-AUC: %.1f%%",
                         bp, info$disease, info$disease_n, info$healthy,
                         auc_val * 100),
    x = 0.52,
    y = 0.44 - (i - 1) * 0.072
  ))
}

# White box behind labels — bottom-right
box_df <- data.frame(xmin = 0.51, xmax = 1.00, ymin = 0.00, ymax = 0.50)

p <- ggplot(roc_data_all, aes(x = FPR, y = TPR, color = BioProject)) +
  
  # Reference diagonal
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey55", linewidth = 0.55) +
  
  # ROC lines
  geom_line(linewidth = 1.5, alpha = 0.92) +
  
  # White AUC table box
  geom_rect(
    data        = box_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "white", color = "#bbbbbb", linewidth = 0.4
  ) +
  
  # AUC labels inside box
  geom_text(
    data     = auc_labels,
    aes(x = x, y = y, label = Full_Label, color = BioProject),
    hjust    = 0,
    size     = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  scale_color_manual(values = colors) +
  
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01),
                     labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0.01),
                     labels = function(x) sprintf("%.1f", x)) +
  
  labs(
    title    = "LOSO Cross-Validation ROC Curves: SVM-RBF-PFI Generalisation Across Six Independent Studies",
    x        = "False Positive Rate",
    y        = "True Positive Rate",
    color    = "Study"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5,
                                     color = "#1a1a1a", margin = margin(b = 5),
                                     lineheight = 1.3),

    axis.title.x      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(t = 8)),
    axis.title.y      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(r = 8)),
    axis.text         = element_text(size = 11, face = "bold", color = "#000000"),
    axis.line         = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks        = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.length = unit(0.12, "cm"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.position   = "none",
    plot.margin       = margin(14, 14, 12, 14)
  )

print(p)

ggsave("2026_Code_Publications/ML_Pictures/loso_roc_curves_six_studies.png",
       plot = p, width = 12, height = 8, dpi = 600, bg = "white")
ggsave("2026_Code_Publications/ML_Pictures/loso_roc_curves_six_studies.pdf",
       plot = p, width = 12, height = 8, device = "pdf", bg = "white")














# Training Dataset: Phenotype Distribution Across Studies - Publication Figure

library(ggplot2)
library(dplyr)
library(forcats)
library(stringr)

bioproject_summary <- read.csv("training_metadata_100K_feb_2026.tsv", sep = "\t")

phenotype_summary <- bioproject_summary %>%
  group_by(Phenotype, Full_Name) %>%
  summarise(Sample_Count = n(), BioProjects = n_distinct(BioProject), .groups = "drop") %>%
  mutate(
    Full_Name = case_when(
      Full_Name == "Atherosclerotic cardiovascular disease" ~ "Atherosclerotic\ncardiovascular disease",
      Full_Name == "Age-related macular degeneration"       ~ "Age-related macular\ndegeneration",
      Full_Name == "Gestational diabetes mellitus"          ~ "Gestational diabetes\nmellitus",
      Full_Name == "Obesity"                                ~ "Obesity",
      TRUE ~ Full_Name
    ),
    Phenotype = case_when(
      Phenotype == "Obesity" ~ "OB",
      TRUE ~ Phenotype
    ),
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
  "Healthy (HC)"                                    = "#4ECDC4",
  "Colorectal Cancer (CRC)"                         = "#E63946",
  "Antibiotic exposer (ABx)"                        = "#1D3557",
  "Parkinson (PD)"                                  = "#9B59B6",
  "Obesity (OB)"                                    = "#2ECC71",
  "Crohn Disease (CD)"                              = "#F39C12",
  "Ulcerative Colitis (UC)"                         = "#E67E22",
  "Atherosclerotic\ncardiovascular disease (ACVD)"  = "#95A5A6",
  "Acute Pancreatitis (AP)"                         = "#16A085",
  "Liver cirrhosis (LV)"                            = "#8E44AD",
  "Ankylosing spondylitis (AS)"                     = "#27AE60",
  "Gastric Cancer (GC)"                             = "#E91E63",
  "Age-related macular\ndegeneration (AMD)"         = "#D68910",
  "Gestational diabetes\nmellitus (GDM)"            = "#A04000"
)

plot_data <- plot_data %>%
  mutate(
    Label_With_Info = sprintf("%s\nn=%d", Phenotype_Label, Sample_Count),
    Label_With_Info = str_wrap(Label_With_Info, width = 35),
    Label_With_Info = fct_reorder(Label_With_Info, Sample_Count, .desc = FALSE)
  )

max_sample_count <- max(plot_data$Sample_Count)
x_limit <- max_sample_count * 1.28

p1 <- ggplot(plot_data, aes(x = Sample_Count, y = Label_With_Info, fill = Phenotype_Label_Display)) +
  geom_bar(stat = "identity", color = "#000000", linewidth = 0.4) +
  geom_text(
    data = plot_data,
    aes(x = Sample_Count, y = Label_With_Info,
        label = sprintf("%d → %.1f%%", BioProjects, Percentage)),
    hjust = -0.1, vjust = 0.5,
    color = "#000000", fontface = "bold", size = 6,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = pheno_colors) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, x_limit),
    breaks = seq(0, max_sample_count, by = 500)
  ) +
  labs(x= "# of Samples",y= "",
    title   = "7,452 Stool-Derived Shotgun Metagenomes Across 13 Phenotypes \nfrom 32 Independent Studies") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 16, color = "#1a1a1a",  # centered
                                     margin = margin(b = 4), lineheight = 1.3),
    axis.text.y       = element_text(size = 14, color = "#000000", face = "bold",hjust = 1, lineheight = 1.1),
    axis.text.x       = element_text(size = 14, color = "#000000", face = "bold"),
    axis.title.x      = element_text(size = 12, face = "bold", color = "#000000",
                                     margin = margin(t = 10)),
    axis.line.x       = element_line(color = "#000000", linewidth = 0.6),
    axis.line.y       = element_blank(),
    axis.ticks.x      = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.y      = element_blank(),
    axis.ticks.length.x = unit(0.15, "cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.border       = element_rect(color = "#000000", fill = NA, linewidth = 0.8),
    plot.margin        = margin(15, 20, 12, 15),
    legend.position    = "none"
  )

print(p1)

# Save as PNG
ggsave("2026_Code_Publications/ML_Pictures/training_dataset_phenotype_distribution.png",
       plot = p1, width = 14, height = 9, dpi = 600, bg = "white")
ggsave("2026_Code_Publications/ML_Pictures/training_dataset_phenotype_distribution.pdf",
       plot = p1, width = 14, height = 9, device = "pdf", bg = "white")



# Validation Dataset: Phenotype Distribution - Publication Figure

library(ggplot2)
library(dplyr)
library(forcats)
library(stringr)

bioproject_summary <- read.csv("fina_validaation_metadata_for_publications.tsv", sep = "\t")

phenotype_summary <- bioproject_summary %>%
  group_by(Phenotype, Full_Name) %>%
  summarise(Sample_Count = n(),
            BioProjects  = n_distinct(BioProject),
            .groups = "drop") %>%
  mutate(
    Full_Name = case_when(
      Full_Name == "Clostridium difficile infection" ~ "Clostridium difficile\ninfection",
      Full_Name == "Type 2 diabetes"                 ~ "Type 2 diabetes",
      Full_Name == "Colorectal Cancer"               ~ "Colorectal Cancer",
      TRUE ~ Full_Name
    ),
    Phenotype = case_when(
      Phenotype == "Obesity" ~ "OB",
      TRUE ~ Phenotype
    ),
    Phenotype_Label = paste0(Full_Name, " (", Phenotype, ")")
  )

plot_data <- phenotype_summary %>%
  select(-Phenotype) %>%
  mutate(
    Phenotype_Label = ifelse(Phenotype_Label == "Healthy (Healthy)",
                             "Healthy (HC)", Phenotype_Label),
    Percentage              = round(Sample_Count / sum(Sample_Count) * 100, 1),
    Phenotype_Label_Display = Phenotype_Label
  )

pheno_colors <- c(
  "Healthy (HC)"                                    = "#4ECDC4",
  "Colorectal Cancer (CRC)"                         = "#E63946",
  "Parkinson (PD)"                                  = "#9B59B6",
  "Obesity (OB)"                                    = "#2ECC71",
  "Crohn Disease (CD)"                              = "#F39C12",
  "Atherosclerotic cardiovascular disease (ACVD)"   = "#95A5A6",
  "Clostridium difficile\ninfection (CDI)"          = "#16A085",
  "Type 2 diabetes (T2D)"                           = "#FF1493"
)

plot_data <- plot_data %>%
  mutate(
    Label_With_Info = sprintf("%s\nn=%d", Phenotype_Label, Sample_Count),
    Label_With_Info = str_wrap(Label_With_Info, width = 35),
    Label_With_Info = fct_reorder(Label_With_Info, Sample_Count, .desc = FALSE)
  )

max_sample_count <- max(plot_data$Sample_Count)
x_limit          <- max_sample_count * 1.32

total_samples    <- sum(plot_data$Sample_Count)
total_bioprojects <- n_distinct(bioproject_summary$BioProject)

p1 <- ggplot(plot_data, aes(x = Sample_Count, y = Label_With_Info,
                            fill = Phenotype_Label_Display)) +
  
  geom_bar(stat = "identity", color = "#000000", linewidth = 0.45) +
  
  geom_text(
    data        = plot_data,
    aes(x     = Sample_Count,
        y     = Label_With_Info,
        label = sprintf("%d → %.1f%%", BioProjects, Percentage)),
    hjust       = -0.08,
    vjust       = 0.5,
    color       = "#000000",
    fontface    = "bold",
    size        = 6,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(values = pheno_colors) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, x_limit),
    breaks = seq(0, max_sample_count, by = 50)
  ) +
  
  labs(
    x     = "# of Samples",
    y     = "",
    title = sprintf("%d Stool-Derived Shotgun Metagenomes Across 7 Phenotypes\nfrom %d Independent Studies",
      total_samples, total_bioprojects
    )
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title          = element_text(hjust = 0.5, face = "bold", size = 16,
                                       color = "#1a1a1a", margin = margin(b = 12),
                                       lineheight = 1.35),
    axis.text.y         = element_text(size = 14, color = "#000000", face = "bold",
                                       hjust = 1, lineheight = 1.15),
    axis.text.x         = element_text(size = 14, color = "#000000", face = "bold"),
    axis.title.x        = element_text(size = 14, face = "bold", color = "#000000",
                                       margin = margin(t = 10)),
    axis.line.x         = element_line(color = "#000000", linewidth = 0.6),
    axis.line.y         = element_blank(),
    axis.ticks.x        = element_line(color = "#000000", linewidth = 0.5),
    axis.ticks.y        = element_blank(),
    axis.ticks.length.x = unit(0.15, "cm"),
    panel.grid.major.y  = element_blank(),
    panel.grid.major.x  = element_line(color = "#eeeeee", linewidth = 0.3),
    panel.grid.minor    = element_blank(),
    panel.background    = element_rect(fill = "white", color = NA),
    plot.background     = element_rect(fill = "white", color = NA),
    panel.border        = element_rect(color = "#000000", fill = NA, linewidth = 0.8),
    plot.margin         = margin(18, 24, 14, 16),
    legend.position     = "none"
  )

print(p1)

# Save as PNG
ggsave("2026_Code_Publications/ML_Pictures/phenotype_distribution_validation_dataset.png",
       plot = p1, width = 15, height = 9, dpi = 600, bg = "white")

# Save as PDF (vector)
ggsave("2026_Code_Publications/ML_Pictures/phenotype_distribution_validation_dataset.pdf",
       plot = p1, width = 15, height = 9, device = "pdf", bg = "white")





########Read and species analysis for the bioprojects.....
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(scales)

df <- read_delim("read_count_training_datasets_with_metadata.csv", delim = "\t")

df <- df %>%
  mutate(
    Total_Reads = Classified + Unclassified,
    Pct_Classified = (Classified / Total_Reads) * 100
  )

bio_avg <- df %>%
  group_by(BioProject) %>%
  summarise(
    avg_classified = mean(Classified, na.rm = TRUE),
    avg_unclassified = mean(Unclassified, na.rm = TRUE),
    avg_total = mean(Classified + Unclassified, na.rm = TRUE),
    total_classified = sum(Classified, na.rm = TRUE),
    total_reads = sum(Classified + Unclassified, na.rm = TRUE),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    BioProject_label = paste0(BioProject, " (n=", count, ")"),
    pct_classified = (total_classified / total_reads) * 100
  ) %>%
  arrange(avg_total)

bio_avg_long <- bio_avg %>%
  select(BioProject_label, avg_classified, avg_unclassified, avg_total, pct_classified) %>%
  pivot_longer(
    cols = c(avg_classified, avg_unclassified),
    names_to = "Category",
    values_to = "avg_reads"
  ) %>%
  mutate(
    Category = factor(Category, 
                      levels = c("avg_unclassified", "avg_classified"),
                      labels = c("Unclassified reads", "Classified reads"))
  )

bio_avg_long$BioProject_label <- factor(bio_avg_long$BioProject_label, 
                                        levels = unique(bio_avg$BioProject_label))

p <- ggplot(bio_avg_long, aes(x = avg_reads / 1e6, y = BioProject_label, fill = Category)) +
  geom_col(data = subset(bio_avg_long, Category == "Classified reads"), 
           position = "identity", width = 0.7, alpha = 0.85) +
  geom_col(data = subset(bio_avg_long, Category == "Unclassified reads"), 
           position = "identity", width = 0.7, alpha = 0.85) +
  geom_text(data = bio_avg,
            aes(x = avg_total / 1e6, y = BioProject_label, 
                label = paste0(round(pct_classified, 1), "%")),
            inherit.aes = FALSE,
            hjust = -0.1,
            color = "black",
            size = 3,
            fontface = "bold") +
  annotate("text", x = max(bio_avg$avg_total / 1e6) * 0.98, 
           y = length(unique(bio_avg$BioProject_label)) + 0.5,
           label = "", size = 3.5, fontface = "bold", hjust = 1) +
  scale_fill_manual(
    name = "Read Type",
    values = c("Unclassified reads" = "#E57373", "Classified reads" = "#66BB6A"),
    breaks = c("Unclassified reads", "Classified reads"),
    labels = c("Unclassified reads", 
               "Classified reads")
  ) +
  scale_x_continuous(
    labels = label_number(suffix = "M", scale = 1),
    expand = expansion(mult = c(0, 0.1)),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Average Classified and Unclassified Reads per BioProject",
    subtitle = "Each bar shows the mean number of reads per sample within each BioProject\n percentages indicate classification rate",
    x = "Average number of reads per sample (millions)",
    y = "BioProject ID (n = number of samples)",
    caption = "Note: Classified reads are taxonomically assigned; Unclassified reads lack taxonomic assignment"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                              margin = margin(b = 5), color = "black"),
    plot.subtitle = element_text(size = 10, hjust = 0.5, 
                                 margin = margin(b = 15), color = "grey30", lineheight = 1.2),
    plot.caption = element_text(size = 9, hjust = 0, color = "grey40", 
                                margin = margin(t = 10), face = "italic"),
    
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
    legend.box.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.7, "cm"),
    legend.margin = margin(8, 10, 8, 10),
    legend.box.margin = margin(0, 0, 10, 0),
    
    plot.margin = margin(15, 30, 15, 15)
  )

print(p)

ggsave("2026_Code_Publications/bioproject_reads_publication.png", 
       plot = p, 
       width = 13, 
       height = 10, 
       dpi = 600, 
       bg = "white")



# Calculate mean reads per BioProject
bio_mean_summary <- df %>%
  group_by(BioProject) %>%
  summarise(
    mean_total_reads = mean(Classified + Unclassified, na.rm = TRUE),
    mean_classified_reads = mean(Classified, na.rm = TRUE),
    mean_unclassified_reads = mean(Unclassified, na.rm = TRUE),
    count = n()
  ) %>%
  mutate(
    pct_classified = (mean_classified_reads / mean_total_reads) * 100,
    BioProject_label = paste0(BioProject, " (", count, ")")
  ) %>%
  arrange(desc(mean_total_reads))
bio_mean_summary




library(dplyr)
library(ggplot2)
library(readr)
library(scales)

df <- read_delim("read_count_training_datasets_with_metadata.csv", delim = "\t")

df <- df %>%
  mutate(
    Total_Reads = Classified + Unclassified,
    Pct_Classified = (Classified / Total_Reads) * 100
  )

cat(strrep("=", 80), "\n")
cat("OVERALL STATISTICS (Across All Samples)\n")
cat(strrep("=", 80), "\n")

avg_classified_per_sample <- mean(df$Classified, na.rm = TRUE)
avg_species_per_sample <- mean(df$Species, na.rm = TRUE)
total_samples <- nrow(df)

cat("Total Number of Samples:", total_samples, "\n")
cat("Average Classified Reads per Sample:", format(round(avg_classified_per_sample), big.mark = ","), "\n")
cat("Average Species Count per Sample:", format(round(avg_species_per_sample, 2), big.mark = ","), "\n\n")

cat(strrep("=", 80), "\n")
cat("STATISTICS BY BIOPROJECT\n")
cat(strrep("=", 80), "\n")

avg_per_bioproject <- df %>%
  group_by(BioProject) %>%
  summarise(
    Avg_Classified_Reads = mean(Classified, na.rm = TRUE),
    Avg_Species_Count = mean(Species, na.rm = TRUE),
    Sample_Count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Avg_Classified_Reads))

avg_per_bioproject <- avg_per_bioproject %>%
  mutate(
    Avg_Classified_Formatted = format(round(Avg_Classified_Reads), big.mark = ","),
    Avg_Species_Formatted = format(round(Avg_Species_Count, 2), big.mark = ",")
  )

cat("\nAverage Statistics per BioProject:\n")
print(avg_per_bioproject %>% select(BioProject, Sample_Count, Avg_Classified_Formatted, Avg_Species_Formatted), n = Inf)

cat("\n")
cat(strrep("=", 80), "\n")
cat("SUMMARY STATISTICS\n")
cat(strrep("=", 80), "\n")
summary_stats <- df %>%
  select(Species, Classified) %>%
  summary()
print(summary_stats)

phenotype_counts <- df %>%
  group_by(BioProject, Phenotype) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(BioProject) %>%
  summarise(
    Phenotype_Details = paste0(Phenotype, " (", Count, ")", collapse = ", "),
    .groups = "drop"
  )

stats_classified <- df %>%
  group_by(BioProject) %>%
  summarise(
    Mean_Classified = mean(Classified, na.rm = TRUE),
    SD_Classified = sd(Classified, na.rm = TRUE),
    Median_Classified = median(Classified, na.rm = TRUE),
    Min_Classified = min(Classified, na.rm = TRUE),
    Max_Classified = max(Classified, na.rm = TRUE),
    Mean_Species = mean(Species, na.rm = TRUE),
    SD_Species = sd(Species, na.rm = TRUE),
    Median_Species = median(Species, na.rm = TRUE),
    Sample_Count = n(),
    .groups = "drop"
  ) %>%
  left_join(phenotype_counts, by = "BioProject") %>%
  arrange(desc(Sample_Count))

cat("\n")
cat(strrep("=", 80), "\n")
cat("DETAILED BIOPROJECT STATISTICS\n")
cat(strrep("=", 80), "\n")
print(stats_classified, n = Inf)

write.table(stats_classified, "BioProject_Classified_Stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nDetailed statistics saved to: BioProject_Classified_Stats.tsv\n")

cat("\n")
cat(strrep("=", 80), "\n")
cat("CORRELATION ANALYSIS\n")
cat(strrep("=", 80), "\n")

pearson_test <- cor.test(df$Species, df$Classified, method = "pearson")
spearman_test <- cor.test(df$Species, df$Classified, method = "spearman")
lm_model <- lm(Classified ~ Species, data = df)
lm_summary <- summary(lm_model)

cat("Pearson Correlation Coefficient (r):", round(pearson_test$estimate, 4), "\n")
cat("Pearson p-value:", format.pval(pearson_test$p.value), "\n\n")

cat("Spearman Correlation Coefficient (ρ):", round(spearman_test$estimate, 4), "\n")
cat("Spearman p-value:", format.pval(spearman_test$p.value), "\n\n")

cat("Linear Regression R-squared:", round(lm_summary$r.squared, 4), "\n")
cat("Adjusted R-squared:", round(lm_summary$adj.r.squared, 4), "\n\n")

cat("Linear Model Summary:\n")
print(lm_summary)

p1 <- ggplot(df, aes(x = Species, y = Classified)) +
  geom_point(color = "#7570b3", size = 2.2, alpha = 0.5, shape = 16) +
  geom_smooth(method = "lm", aes(color = "Linear (Pearson)"), 
              se = TRUE, fill = "#1b9e77", alpha = 0.15, linewidth = 1.1) +
  geom_smooth(method = "loess", aes(color = "Non-linear (Spearman)"), 
              se = FALSE, linewidth = 1.1, linetype = "dashed") +
  scale_color_manual(
    name = "Regression Type",
    values = c("Linear (Pearson)" = "#1b9e77", "Non-linear (Spearman)" = "#d95f02"),
    labels = c(
      sprintf("Linear (Pearson r = %.3f***)", pearson_test$estimate),
      sprintf("Non-linear (Spearman ρ = %.3f***)", spearman_test$estimate)
    )
  ) +
  scale_y_continuous(labels = label_number(suffix = "M", scale = 1e-6, accuracy = 0.1)) +
  scale_x_continuous(labels = label_comma()) +
  labs(
    x = "Number of Species per Sample",
    y = "Classified Reads (Millions)",
    title = "Correlation Between Species Diversity and Classified Reads",
    subtitle = sprintf("R² = %.3f, n = %d samples", lm_summary$r.squared, nrow(df))
  ) +
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

print(p1)

ggsave("2026_Code_Publications/Species_vs_Classified_Reads.png", 
       plot = p1, 
       width = 10, 
       height = 8, 
       dpi = 600, 
       bg = "white")






