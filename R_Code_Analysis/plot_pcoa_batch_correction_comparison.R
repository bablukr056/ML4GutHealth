# PCoA visualization: batch correction comparison (training and validation, before vs after)

setwd("/Users/bablu/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/CODE")

library(ggplot2)
library(patchwork)
library(readr)
library(dplyr)

output_path <- "pco_result_all"

# Load PCoA coordinates and variance-explained summary
train_before <- read_csv(file.path(output_path, "PCoA_training_before_correction_with_metadata.csv"), show_col_types = FALSE)
train_after  <- read_csv(file.path(output_path, "PCoA_training_after_correction_with_metadata.csv"),  show_col_types = FALSE)
val_before   <- read_csv(file.path(output_path, "PCoA_validation_before_correction_with_metadata.csv"), show_col_types = FALSE)
val_after    <- read_csv(file.path(output_path, "PCoA_validation_after_correction_with_metadata.csv"),  show_col_types = FALSE)
var_exp      <- read_csv(file.path(output_path, "PCoA_variance_explained_summary.csv"), show_col_types = FALSE)

make_pcoa_plot <- function(df, var_row, title) {
  pco1_lab <- sprintf("PCo1 (%.2f%%)", var_row$PCo1_variance_explained * 100)
  pco2_lab <- sprintf("PCo2 (%.2f%%)", var_row$PCo2_variance_explained * 100)

  ggplot(df, aes(x = PCo1, y = PCo2, color = BioProject)) +
    geom_point(size = 1.4, alpha = 0.7, show.legend = FALSE) +
    labs(title = title, x = pco1_lab, y = pco2_lab) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0, margin = margin(b = 8)),
      axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 8)),
      axis.text = element_text(size = 11, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
      plot.margin = margin(10, 10, 10, 10)
    )
}

p1 <- make_pcoa_plot(
  train_before, var_exp[var_exp$Dataset == "Training_Before", ],
  "Training Set \u2014 Before Batch Correction"
)
p2 <- make_pcoa_plot(
  train_after, var_exp[var_exp$Dataset == "Training_After", ],
  "Training Set \u2014 After Batch Correction"
)
p3 <- make_pcoa_plot(
  val_before, var_exp[var_exp$Dataset == "Validation_Before", ],
  "Validation Set \u2014 Before Batch Correction"
)
p4 <- make_pcoa_plot(
  val_after, var_exp[var_exp$Dataset == "Validation_After", ],
  "Validation Set \u2014 After Batch Correction"
)

final_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

final_plot

# Export figure (PNG, PDF, TIFF)
ggsave(
  file.path(output_path, "PCoA_batch_correction_comparison_R.png"),
  plot = final_plot, width = 14, height = 12, dpi = 600, units = "in"
)

ggsave(
  file.path(output_path, "PCoA_batch_correction_comparison_R.pdf"),
  plot = final_plot, width = 14, height = 12, units = "in"
)

ggsave(
  file.path(output_path, "PCoA_batch_correction_comparison_R.tiff"),
  plot = final_plot, width = 14, height = 12, dpi = 600, units = "in",
  compression = "lzw"
)
