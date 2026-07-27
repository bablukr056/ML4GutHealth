setwd("/Users/bablu/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/CODE/Thesis_Figures_R")

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)


summary_df <- read_csv("Summary_Table_WITHtuning_all_configurations.csv", show_col_types = FALSE)
summary_df$Configuration <- factor(summary_df$Configuration, levels = summary_df$Configuration)

plot_df <- summary_df %>%
  select(Configuration, F1, ROC_AUC) %>%
  pivot_longer(cols = c(F1, ROC_AUC), names_to = "Metric", values_to = "Score") %>%
  mutate(Metric = recode(Metric, "F1" = "Test Set F1 Score", "ROC_AUC" = "Test Set ROC-AUC"))

config_colors <- c("#3182bd", "#e6550d", "#31a354", "#756bb1", "#636363")
names(config_colors) <- levels(summary_df$Configuration)

theme_publication <- theme_bw(base_size = 15) +
  theme(
    plot.title       = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title.x     = element_text(face = "bold", size = 15, margin = margin(t = 10)),
    axis.title.y     = element_text(face = "bold", size = 15, margin = margin(r = 10)),
    axis.text.x      = element_text(face = "bold", size = 12, angle = 30, hjust = 1, color = "black"),
    axis.text.y      = element_text(face = "bold", size = 12, color = "black"),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed", linewidth = 0.3, color = "grey70"),
    panel.border     = element_rect(color = "black", linewidth = 0.8),
    strip.text       = element_text(face = "bold", size = 13),
    plot.margin      = margin(10, 10, 10, 10)
  )

make_metric_plot <- function(metric_label) {
  df_sub <- plot_df %>% filter(Metric == metric_label)
  
  ggplot(df_sub, aes(x = Configuration, y = Score, fill = Configuration)) +
    geom_col(color = "black", linewidth = 0.6, width = 0.65, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.3f", Score)),
              vjust = -0.5, fontface = "bold", size = 4.2) +
    scale_fill_manual(values = config_colors) +
    scale_y_continuous(limits = c(0, 1.08), expand = expansion(mult = c(0, 0.02))) +
    labs(title = paste0(metric_label, " — All Configurations"),
         x = "Configuration", y = metric_label) +
    theme_publication
}

p_f1  <- make_metric_plot("Test Set F1 Score")
p_auc <- make_metric_plot("Test Set ROC-AUC")

combined_plot <- p_f1 + p_auc + plot_layout(ncol = 2)
combined_plot

out_base <- "Figure_WITHtuning_ALL_configurations_combined_R"

ggsave(paste0(out_base, ".png"), combined_plot, width = 13, height = 6, dpi = 600, bg = "white")
ggsave(paste0(out_base, ".pdf"), combined_plot, width = 13, height = 6, device = cairo_pdf, bg = "white")
ggsave(paste0(out_base, ".tiff"), combined_plot, width = 13, height = 6, dpi = 600,
       compression = "lzw", bg = "white")
