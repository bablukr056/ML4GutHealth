# ==============================================================================
# Beta Diversity Visualization (PCoA, Aitchison Distance)
#
# This script generates the beta-diversity figures used in the manuscript.
# PCoA ordination was performed in Python (scikit-bio) on Aitchison distances
# (Euclidean distance of CLR-transformed species-level abundance data); the
# resulting PCoA coordinates were exported and are visualized here in R.
#
# Three figures are produced:
#   1. PCoA by health status (healthy vs. non-healthy)
#   2. PCoA by disease phenotype (six selected conditions)
#   3. A combined two-panel figure (health status + phenotype)
# ==============================================================================

library(ggplot2)
library(ggforce)
library(dplyr)
library(scales)
library(patchwork)

data <- read.delim("beta_diversity_pcoa_scores_all_beta.tsv", sep = "\t")
cat("Data dimensions:", dim(data)[1], "samples x", dim(data)[2], "columns\n")

pc1_col <- grep("^PC1", names(data), value = TRUE)[1]
pc2_col <- grep("^PC2", names(data), value = TRUE)[1]

if (is.na(pc1_col) || is.na(pc2_col)) {
  stop("Could not find PC1/PC2 columns automatically. Check colnames(data) and set pc1_col/pc2_col manually.")
}
cat(sprintf("Using columns: %s (x-axis), %s (y-axis)\n", pc1_col, pc2_col))

dir.create("beta_diversity_figure", showWarnings = FALSE)

theme_manuscript <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.text            = element_text(size = 11, color = "black", face = "bold"),
      axis.title.x         = element_text(size = 13, color = "black", face = "bold", margin = margin(t = 8, b = 4)),
      axis.title.y         = element_text(size = 13, color = "black", face = "bold", margin = margin(r = 8, l = 4)),
      plot.title           = element_text(size = 14.5, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 10)),
      panel.grid.major.y   = element_line(color = "#E5E5E5", linewidth = 0.3),
      panel.grid.minor.y   = element_line(color = "#F5F5F5", linewidth = 0.2),
      panel.grid.major.x   = element_blank(),
      panel.grid.minor.x   = element_blank(),
      plot.margin          = margin(14, 16, 12, 12),
      panel.background     = element_rect(fill = "white", color = NA),
      plot.background      = element_rect(fill = "white", color = NA),
      panel.border         = element_rect(color = "#333333", fill = NA, linewidth = 0.7),
      legend.position      = "top",
      legend.title         = element_text(size = 11.5, face = "bold", color = "black"),
      legend.text          = element_text(size = 10.5, face = "bold", color = "black"),
      legend.background    = element_rect(fill = "white", color = NA),
      legend.key           = element_rect(fill = "white", color = NA),
      legend.key.size      = unit(0.8, "lines"),
      axis.line            = element_line(color = "#333333", linewidth = 0.5),
      axis.ticks           = element_line(color = "#333333", linewidth = 0.5),
      axis.ticks.length    = unit(0.15, "cm")
    )
}


# ==============================================================================
# Figure 1: PCoA by health status
# ==============================================================================

health_colors <- c("HEALTHY" = "#56B4E9", "NON-HEALTHY" = "#fc8d62")
health_data <- data[data$Health_status %in% names(health_colors), ]

health_counts <- table(health_data$Health_status)
health_labels <- setNames(paste0(names(health_counts), " (n=", health_counts, ")"), names(health_counts))

hx_range <- range(health_data[[pc1_col]], na.rm = TRUE)
hy_range <- range(health_data[[pc2_col]], na.rm = TRUE)
h_label_x <- hx_range[2] - 0.02 * diff(hx_range)
h_label_y <- hy_range[2] + 0.06 * diff(hy_range)

p_health <- ggplot(health_data, aes(x = .data[[pc1_col]], y = .data[[pc2_col]], color = Health_status)) +
  geom_point(size = 1, alpha = 0.5, shape = 16) +
  stat_ellipse(aes(group = Health_status), type = "norm", level = 0.95, linetype = "dashed", linewidth = 1) +
  stat_summary(aes(group = Health_status), fun = "mean", geom = "point", shape = 16, size = 5.5, color = "black") +
  stat_summary(aes(group = Health_status), fun = "mean", geom = "point", shape = 16, size = 4.5) +
  scale_color_manual(values = health_colors, labels = health_labels) +
  annotate("text", x = h_label_x, y = h_label_y,
           label = "PERMANOVA: pseudo-F = 41.89, p = 0.001, R\u00b2 = 0.0056",
           hjust = 1, vjust = 1, size = 4.2, fontface = "bold") +
  labs(x = "PC1 (8.73%)", y = "PC2 (7.72%)",
       title = "PCoA of Gut Microbial Composition (Aitchison Distance)",
       color = "Health Status") +
  theme_manuscript()

ggsave(file.path("beta_diversity_figure", "pcoa_healthy_vs_nonhealthy.png"),
       plot = p_health, width = 11, height = 8.5, dpi = 600, bg = "white")
ggsave(file.path("beta_diversity_figure", "pcoa_healthy_vs_nonhealthy.pdf"),
       plot = p_health, width = 11, height = 8.5, bg = "white")
ggsave(file.path("beta_diversity_figure", "pcoa_healthy_vs_nonhealthy.tiff"),
       plot = p_health, width = 11, height = 8.5, dpi = 600, compression = "lzw", bg = "white")


# ==============================================================================
# Figure 2: PCoA by disease phenotype (six selected conditions)
# ==============================================================================

selected_conditions <- c("Healthy", "CD", "LV", "CRC", "PD", "Obesity")

pheno_colors <- c(
  "Healthy" = "#56B4E9",
  "CRC"     = "#D55E00",
  "PD"      = "#CC79A7",
  "Obesity" = "#009E73",
  "CD"      = "#F0E442",
  "LV"      = "#7570B3"
)

pheno_data <- data %>%
  dplyr::mutate(Phenotype = trimws(Phenotype)) %>%
  dplyr::filter(Phenotype %in% selected_conditions) %>%
  dplyr::mutate(Phenotype = factor(Phenotype, levels = names(pheno_colors)))

pheno_counts <- table(pheno_data$Phenotype)
pheno_labels <- setNames(paste0(names(pheno_colors), " (n=", pheno_counts[names(pheno_colors)], ")"), names(pheno_colors))

px_range <- range(pheno_data[[pc1_col]], na.rm = TRUE)
py_range <- range(pheno_data[[pc2_col]], na.rm = TRUE)
p_label_x <- px_range[2] - 0.02 * diff(px_range)
p_label_y <- py_range[2] + 0.06 * diff(py_range)

p_pheno <- ggplot(pheno_data, aes(x = .data[[pc1_col]], y = .data[[pc2_col]], color = Phenotype)) +
  geom_point(size = 1, alpha = 0.5, shape = 16) +
  stat_ellipse(aes(group = Phenotype), type = "norm", level = 0.95, linetype = "dashed", linewidth = 1) +
  stat_summary(aes(group = Phenotype), fun = "mean", geom = "point", shape = 16, size = 5.5, color = "black") +
  stat_summary(aes(group = Phenotype), fun = "mean", geom = "point", shape = 16, size = 4.5) +
  scale_color_manual(values = pheno_colors, labels = pheno_labels) +
  annotate("text", x = p_label_x, y = p_label_y,
           label = "PERMANOVA: pseudo-F = 26.09, p = 0.001, R\u00b2 = 0.0404",
           hjust = 1, vjust = 1, size = 4.2, fontface = "bold") +
  labs(x = "PC1 (8.73%)", y = "PC2 (7.72%)",
       title = "Gut Microbiome Beta Diversity by Phenotype",
       color = "Phenotype") +
  theme_manuscript()

ggsave(file.path("beta_diversity_figure", "pcoa_selected_conditions.png"),
       plot = p_pheno, width = 11, height = 8.5, dpi = 600, bg = "white")
ggsave(file.path("beta_diversity_figure", "pcoa_selected_conditions.pdf"),
       plot = p_pheno, width = 11, height = 8.5, bg = "white")
ggsave(file.path("beta_diversity_figure", "pcoa_selected_conditions.tiff"),
       plot = p_pheno, width = 11, height = 8.5, dpi = 600, compression = "lzw", bg = "white")


# ==============================================================================
# Figure 3: Combined two-panel figure (health status + phenotype)
# ==============================================================================

combined_fig <- p_health / p_pheno +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

print(combined_fig)

ggsave(file.path("beta_diversity_figure", "pcoa_combined_panels.png"),
       plot = combined_fig, width = 11, height = 16, dpi = 600, bg = "white")
ggsave(file.path("beta_diversity_figure", "pcoa_combined_panels.pdf"),
       plot = combined_fig, width = 11, height = 16, bg = "white")
ggsave(file.path("beta_diversity_figure", "pcoa_combined_panels.tiff"),
       plot = combined_fig, width = 11, height = 16, dpi = 600, compression = "lzw", bg = "white")

cat("\nAll beta-diversity figures saved to: beta_diversity_figure/\n")
