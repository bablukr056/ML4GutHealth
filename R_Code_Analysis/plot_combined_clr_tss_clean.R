#!/usr/bin/env Rscript

# ==============================================================================
# rerun_normalizations_summary.csv 
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

summary_file <- "/Users/bablu/OneDrive - Università degli Studi di Milano/CODE/code_submission_research_paper/rerun_normalizations_summary.csv"
output_dir <- dirname(summary_file)

cat("Reading summary table from:", summary_file, "\n")
df <- read.csv(summary_file, stringsAsFactors = FALSE)
print(df)

# ------------------------------------------------------------------------
# Readable display labels for each configuration
# ------------------------------------------------------------------------
display_label_map <- c(
  "Raw_Counts" = "Raw counts",
  "TSS_Relative_Abundance" = "TSS\n(relative abundance)",
  "CLR_pseudocount_1.0" = "CLR\npseudocount = 1.0\n(manuscript)",
  "CLR_pseudocount_0.5" = "CLR\npseudocount = 0.5",
  "CLR_pseudocount_0.1" = "CLR\npseudocount = 0.1",
  "CLR_multiplicative_replacement" = "CLR\nmultiplicative\nreplacement"
)

df$Display_Label <- display_label_map[df$Configuration]
df$Display_Label[is.na(df$Display_Label)] <- df$Configuration[is.na(df$Display_Label)]

df$Display_Label <- factor(df$Display_Label, levels = df$Display_Label)
df$Highlight <- ifelse(df$Configuration == "CLR_pseudocount_1.0", "Manuscript configuration", "Comparison")

palette_colors <- c("#969696", "#fdae6b", "#3182bd", "#6baed6", "#9ecae1", "#c994c7")
palette_colors <- palette_colors[seq_len(nrow(df))]
names(palette_colors) <- df$Display_Label

theme_annotated <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(size = 9.5, color = "black", face = "bold",
                                 angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10.5, color = "black", face = "bold"),
      axis.title = element_text(face = "bold", size = 12.5),
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position = "none",
      plot.margin = margin(10, 12, 15, 10)
    )
}

make_metric_plot <- function(data, metric_col, metric_label) {
  data$Value <- data[[metric_col]]
  data$LineWidth <- ifelse(data$Highlight == "Manuscript configuration", 1.6, 0.6)
  data$LineColor <- ifelse(data$Highlight == "Manuscript configuration", "red", "black")
  
  ggplot(data, aes(x = Display_Label, y = Value, fill = Display_Label)) +
    geom_col(color = data$LineColor, linewidth = data$LineWidth, width = 0.65) +
    geom_text(aes(label = sprintf("%.1f%%", Value * 100)),
              vjust = -0.6, fontface = "bold", size = 3.8) +
    scale_fill_manual(values = palette_colors) +
    scale_y_continuous(
      limits = c(0, 1.08),
      expand = c(0, 0),
      labels = label_percent(accuracy = 1)
    ) +
    labs(
      title = paste("Test Set", metric_label, "by Normalization / Zero-Handling Strategy"),
      x = NULL,
      y = metric_label
    ) +
    theme_annotated()
}

p_f1 <- make_metric_plot(df, "F1", "F1 Score")
p_auc <- make_metric_plot(df, "ROC_AUC", "ROC-AUC")

final_plot <- (p_f1 | p_auc) +
  plot_annotation(
    caption = expression(bold("Red outline") ~ "denotes the normalization strategy (CLR, pseudocount = 1.0) used in the primary manuscript analysis."),
    theme = theme(plot.caption = element_text(size = 10, hjust = 0.5))
  )

final_plot

ggsave(file.path(output_dir, "Figure_normalization_zerohandling_sensitivity_comparison.png"),
       plot = final_plot, width = 14, height = 6, dpi = 600, units = "in", bg = "white")
ggsave(file.path(output_dir, "Figure_normalization_zerohandling_sensitivity_comparison.pdf"),
       plot = final_plot, width = 14, height = 6, units = "in", bg = "white")
ggsave(file.path(output_dir, "Figure_normalization_zerohandling_sensitivity_comparison.tiff"),
       plot = final_plot, width = 14, height = 6, dpi = 600, units = "in",
       compression = "lzw", bg = "white")

cat("\nAnnotated combined figure saved to:", output_dir, "\n")
cat("DONE.\n")
