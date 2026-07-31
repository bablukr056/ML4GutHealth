# ============================================================
# Beta diversity PCoA + Alpha diversity analysis
# Figures: PCoA by health status/phenotype; alpha diversity violin
# plots by health status/phenotype (Shannon, Chao1); descriptive
# summary tables of alpha diversity
# ============================================================

setwd("/Users/bablu/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/CODE/Thesis_Figures_R")

# ============================================================
# PART 1: Beta diversity — PCoA (Aitchison distance)
# Panel A: samples colored by Health Status (Healthy vs Non-Healthy)
# Panel B: samples colored by Phenotype (disease category)
# ============================================================

library(ggplot2)
library(ggforce)
library(dplyr)
library(scales)
library(patchwork)

data <- read.delim("beta_diversity_pcoa_scores_all_beta.tsv", sep = "\t")
cat("Data Dimensions:\n"); print(dim(data))
cat("Column Names:\n"); print(colnames(data))

pc1_col <- grep("^PC1", names(data), value = TRUE)[1]
pc2_col <- grep("^PC2", names(data), value = TRUE)[1]

if (is.na(pc1_col) || is.na(pc2_col)) {
  stop("Could not find PC1/PC2 columns automatically. Check colnames(data) above and set pc1_col/pc2_col manually.")
}
cat(sprintf("Using columns: %s (x-axis), %s (y-axis)\n", pc1_col, pc2_col))

theme_thesis <- function(base_size = 13) {
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
      legend.box.margin    = margin(0, 0, 0, 0),
      axis.line            = element_line(color = "#333333", linewidth = 0.5),
      axis.ticks           = element_line(color = "#333333", linewidth = 0.5),
      axis.ticks.length    = unit(0.15, "cm")
    )
}

# --- Panel A: Healthy vs Non-Healthy ---
health_colors <- c("HEALTHY" = "#56B4E9", "NON-HEALTHY" = "#fc8d62")

health_data <- data[data$Health_status %in% names(health_colors), ]

health_counts <- table(health_data$Health_status)
health_labels <- setNames(
  paste0(names(health_counts), " (n=", health_counts, ")"),
  names(health_counts)
)

hx_range <- range(health_data[[pc1_col]], na.rm = TRUE)
hy_range <- range(health_data[[pc2_col]], na.rm = TRUE)
h_label_x <- hx_range[2] - 0.02 * diff(hx_range)
h_label_y <- hy_range[2] + 0.06 * diff(hy_range)

p_health <- ggplot(health_data, aes(x = .data[[pc1_col]], y = .data[[pc2_col]], color = Health_status, fill = Health_status)) +
  geom_point(size = 0.4, alpha = 0.35, shape = 16) +
  geom_mark_ellipse(aes(group = Health_status), expand = unit(2, "mm"),
                    linewidth = 0.9, linetype = "dashed", alpha = 0.08, show.legend = FALSE) +
  stat_summary(aes(group = Health_status), fun = "mean", geom = "point", shape = 16, size = 3.5, color = "black") +
  stat_summary(aes(group = Health_status), fun = "mean", geom = "point", shape = 16, size = 2.7) +
  scale_color_manual(values = health_colors, labels = health_labels) +
  scale_fill_manual(values = health_colors, guide = "none") +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 4, alpha = 1))) +
  annotate("text", x = h_label_x, y = h_label_y,
           label = "PERMANOVA: pseudo-F = 41.89, p = 0.001, R\u00b2 = 0.0056",
           hjust = 1, vjust = 1, size = 4.2, fontface = "bold") +
  labs(x = "PC1 (8.73%)", y = "PC2 (7.72%)",
       title = "PCoA of Gut Microbial Composition (Aitchison Distance)",
       color = "Health Status") +
  theme_thesis()
p_health

# --- Panel B: Phenotype ---
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
pheno_labels <- setNames(
  paste0(names(pheno_counts), " (n=", pheno_counts, ")"),
  names(pheno_counts)
)

px_range <- range(pheno_data[[pc1_col]], na.rm = TRUE)
py_range <- range(pheno_data[[pc2_col]], na.rm = TRUE)
p_label_x <- px_range[2] - 0.02 * diff(px_range)
p_label_y <- py_range[2] + 0.06 * diff(py_range)

p_pheno <- ggplot(pheno_data, aes(x = .data[[pc1_col]], y = .data[[pc2_col]], color = Phenotype, fill = Phenotype)) +
  geom_point(size = 0.4, alpha = 0.35, shape = 16) +
  geom_mark_ellipse(aes(group = Phenotype), expand = unit(2, "mm"),
                    linewidth = 0.9, linetype = "dashed", alpha = 0.08, show.legend = FALSE) +
  stat_summary(aes(group = Phenotype), fun = "mean", geom = "point", shape = 16, size = 3.5, color = "black") +
  stat_summary(aes(group = Phenotype), fun = "mean", geom = "point", shape = 16, size = 2.7) +
  scale_color_manual(values = pheno_colors, labels = pheno_labels) +
  scale_fill_manual(values = pheno_colors, guide = "none") +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 4, alpha = 1))) +
  annotate("text", x = p_label_x, y = p_label_y,
           label = "PERMANOVA: pseudo-F = 26.09, p = 0.001, R\u00b2 = 0.0404",
           hjust = 1, vjust = 1, size = 4.2, fontface = "bold") +
  labs(x = "PC1 (8.73%)", y = "PC2 (7.72%)",
       title = "Gut Microbiome Beta Diversity by Phenotype",
       color = "Phenotype") +
  theme_thesis()
p_pheno

# --- Combine Panel A + Panel B and export ---
combined_fig <- p_health / p_pheno +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

print(combined_fig)

output_dir <- "beta_diversity_figure_submission"
dir.create(output_dir, showWarnings = FALSE)

ggsave(file.path(output_dir, "beta_diversity_pcoa_combined.png"),
       plot = combined_fig, width = 11, height = 16, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "beta_diversity_pcoa_combined.pdf"),
       plot = combined_fig, width = 11, height = 16, device = "pdf", bg = "white")
ggsave(file.path(output_dir, "beta_diversity_pcoa_combined.tiff"),
       plot = combined_fig, width = 11, height = 16, dpi = 600, compression = "lzw", bg = "white")


# ============================================================
# PART 2: Alpha diversity — Shannon and Chao1 by Health Status and Phenotype
# 4-panel figure: (A) Chao1 by Health Status, (B) Chao1 by Phenotype,
#                 (C) Shannon by Health Status, (D) Shannon by Phenotype
# Statistical summary (Wilcoxon test, Cliff's delta) exported as a table
# ============================================================

library(ggplot2)
library(dplyr)
library(ggpubr)
library(effsize)
library(scales)
library(patchwork)

div_df <- read.delim("alpha_diversity_with_metadata_greater100K.csv", sep = "\t", header = TRUE)
dim(div_df)

div_df <- div_df %>%
  mutate(
    BioProject      = as.factor(BioProject),
    Health_status   = factor(Health_status, levels = c("HEALTHY", "NON-HEALTHY")),
    Shannon.Entropy = as.numeric(Shannon.Entropy),
    Chao1           = as.numeric(Chao1),
    Phenotype       = as.factor(Phenotype)
  )

y_labs <- c(Shannon.Entropy = "Shannon Diversity Index",
            Chao1           = "Chao1 Diversity Index")

pheno_colors <- c(
  "Healthy" = "#2E86AB",
  "CRC"     = "#E4572E",
  "PD"      = "#CC79A7",
  "Obesity" = "#009E73",
  "CD"      = "#F0E442",
  "UC"      = "#E69F00",
  "ACVD"    = "#999999",
  "AP"      = "#1B9E77",
  "LV"      = "#7570B3",
  "AS"      = "#66A61E",
  "GC"      = "#E7298A",
  "AMD"     = "#A6761D",
  "GDM"     = "#8C510A"
)

fill_colors    <- c(HEALTHY = "#2E86AB", `NON-HEALTHY` = "#E4572E")
outline_colors <- c(HEALTHY = "#1B5876", `NON-HEALTHY` = "#A8391A")

theme_thesis <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position     = "none",
      axis.text           = element_text(size = 11, color = "grey15", face = "bold"),
      axis.title          = element_text(size = 12.5, color = "black", face = "bold"),
      axis.title.y        = element_text(margin = margin(r = 12)),
      axis.title.x        = element_text(margin = margin(t = 12)),
      axis.ticks          = element_line(color = "grey40", linewidth = 0.4),
      axis.ticks.length   = unit(4, "pt"),
      plot.title          = element_text(size = 13.5, face = "bold", hjust = 0,
                                         margin = margin(b = 10), color = "black"),
      plot.caption        = element_text(size = 8.5, color = "grey45", hjust = 1,
                                         margin = margin(t = 8)),
      panel.grid          = element_blank(),
      panel.background    = element_rect(fill = "white", color = NA),
      plot.background     = element_rect(fill = "white", color = NA),
      panel.border        = element_rect(color = "grey20", fill = NA, linewidth = 0.6),
      plot.margin         = margin(14, 16, 12, 12)
    )
}

# Violin + boxplot of a diversity measure across Phenotype, with Wilcoxon test vs Healthy
create_phenotype_plot <- function(data, measure, y_label) {
  phenotype_counts <- table(data$Phenotype)
  labels <- paste0(names(phenotype_counts), "\n(n = ", phenotype_counts, ")")
  color_palette <- pheno_colors[names(pheno_colors) %in% names(phenotype_counts)]

  ggplot(data, aes(x = Phenotype, y = .data[[measure]], fill = Phenotype)) +
    geom_violin(trim = FALSE, alpha = 0.75, linewidth = 0.5, color = "grey30") +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA,
                 linewidth = 0.4, color = "grey20") +
    stat_compare_means(
      method = "wilcox.test",
      ref.group = "Healthy",
      label = "p.signif",
      p.adjust.method = "fdr",
      label.y = max(data[[measure]], na.rm = TRUE) * 1.08,
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")
      ),
      size = 4.2,
      vjust = -0.5
    ) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(labels = number_format(accuracy = 0.1),
                       expand = expansion(mult = c(0.05, 0.18))) +
    scale_fill_manual(values = color_palette) +
    labs(
      x = "Phenotype",
      y = y_label,
      title = y_label,
      caption = "Wilcoxon rank-sum test vs Healthy, FDR-adjusted; *** p<0.001, ** p<0.01, * p<0.05, ns = not significant"
    ) +
    theme_thesis()
}

# Violin + boxplot of a diversity measure across Health Status, with Wilcoxon test and Cliff's delta
create_health_status_plot <- function(data, measure) {
  f <- reformulate("Health_status", measure)
  wilcox_res <- wilcox.test(f, data = data, conf.int = TRUE)
  cliff_res  <- cliff.delta(f, data = data)
  n  <- table(data$Health_status)
  ym <- max(data[[measure]], na.rm = TRUE)

  plot <- ggplot(data, aes(Health_status, .data[[measure]], fill = Health_status)) +
    geom_violin(aes(color = Health_status), trim = FALSE, alpha = 0.75, linewidth = 0.5) +
    geom_boxplot(width = 0.14, fill = "white", outlier.shape = NA,
                 linewidth = 0.45, color = "grey20") +
    annotate("text", x = 1.5, y = ym * 1.08,
             label = sprintf("P = %s", format.pval(wilcox_res$p.value, digits = 3)),
             size = 3.8, fontface = "bold", color = "grey10") +
    annotate("text", x = 1.5, y = ym * 1.015,
             label = sprintf("Cliff's \u0394 = %.2f (%.2f, %.2f)",
                             cliff_res$estimate, cliff_res$conf.int[1], cliff_res$conf.int[2]),
             size = 3.4, fontface = "italic", color = "grey25") +
    scale_x_discrete(labels = function(x) paste0(x, "\n(n = ", n[x], ")")) +
    scale_y_continuous(labels = number_format(accuracy = 0.1),
                       expand = expansion(mult = c(0.05, 0.18))) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = outline_colors) +
    labs(x = "Health Status", y = y_labs[[measure]],
         title = y_labs[[measure]],
         caption = "Wilcoxon rank-sum test; Cliff's delta effect size") +
    theme_thesis()

  list(plot = plot, wilcox = wilcox_res, cliff = cliff_res, n = n)
}

alpha_data <- div_df[, c("Health_status", "Shannon.Entropy", "Chao1", "Phenotype")]

shannon_pheno  <- create_phenotype_plot(alpha_data, "Shannon.Entropy", y_labs[["Shannon.Entropy"]])
chao1_pheno    <- create_phenotype_plot(alpha_data, "Chao1", y_labs[["Chao1"]])
shannon_health <- create_health_status_plot(alpha_data, "Shannon.Entropy")
chao1_health   <- create_health_status_plot(alpha_data, "Chao1")

combined_fig <- (chao1_health$plot + (chao1_pheno + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))) /
  (shannon_health$plot + (shannon_pheno + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 15, face = "bold"))

print(combined_fig)

dir.create("diversity_figure", showWarnings = FALSE)

ggsave(
  file.path("diversity_figure", "alpha_diversity_combined_panels.png"),
  combined_fig, width = 20, height = 14, dpi = 600, units = "in", bg = "white"
)

ggsave(
  file.path("diversity_figure", "alpha_diversity_combined_panels.pdf"),
  combined_fig, width = 20, height = 14, device = "pdf", bg = "white"
)

# Statistical summary table (Wilcoxon test + Cliff's delta) for Shannon and Chao1
stats_tbl <- bind_rows(
  data.frame(
    Measure        = y_labs[["Shannon.Entropy"]],
    N_Healthy      = shannon_health$n[["HEALTHY"]],
    N_NonHealthy   = shannon_health$n[["NON-HEALTHY"]],
    W              = shannon_health$wilcox$statistic,
    P_value        = shannon_health$wilcox$p.value,
    CI_low         = shannon_health$wilcox$conf.int[1],
    CI_high        = shannon_health$wilcox$conf.int[2],
    Cliffs_delta   = shannon_health$cliff$estimate,
    Cliffs_CI_low  = shannon_health$cliff$conf.int[1],
    Cliffs_CI_high = shannon_health$cliff$conf.int[2]
  ),
  data.frame(
    Measure        = y_labs[["Chao1"]],
    N_Healthy      = chao1_health$n[["HEALTHY"]],
    N_NonHealthy   = chao1_health$n[["NON-HEALTHY"]],
    W              = chao1_health$wilcox$statistic,
    P_value        = chao1_health$wilcox$p.value,
    CI_low         = chao1_health$wilcox$conf.int[1],
    CI_high        = chao1_health$wilcox$conf.int[2],
    Cliffs_delta   = chao1_health$cliff$estimate,
    Cliffs_CI_low  = chao1_health$cliff$conf.int[1],
    Cliffs_CI_high = chao1_health$cliff$conf.int[2]
  )
)

write.table(stats_tbl, file.path("diversity_figure", "alpha_diversity_stats_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# ============================================================
# PART 3: Alpha diversity descriptive summary
# Mean/SD/median of Shannon and Chao1, grouped by Phenotype and by Health Status
# Exported as supplementary tables (not a figure)
# ============================================================

library(dplyr)

div_df <- read.delim("alpha_diversity_with_metadata_greater100K.csv",
                     sep = "\t", stringsAsFactors = FALSE)

descriptive_by_phenotype <- div_df %>%
  group_by(Phenotype) %>%
  summarise(
    N = n(),
    Shannon_mean = round(mean(Shannon.Entropy, na.rm = TRUE), 3),
    Shannon_sd = round(sd(Shannon.Entropy, na.rm = TRUE), 3),
    Shannon_median = round(median(Shannon.Entropy, na.rm = TRUE), 3),
    Chao1_mean = round(mean(Chao1, na.rm = TRUE), 1),
    Chao1_sd = round(sd(Chao1, na.rm = TRUE), 1),
    Chao1_median = round(median(Chao1, na.rm = TRUE), 1)
  ) %>%
  arrange(desc(N))

print(descriptive_by_phenotype)

descriptive_by_health <- div_df %>%
  group_by(Health_status) %>%
  summarise(
    N = n(),
    Shannon_mean = round(mean(Shannon.Entropy, na.rm = TRUE), 3),
    Shannon_sd = round(sd(Shannon.Entropy, na.rm = TRUE), 3),
    Shannon_median = round(median(Shannon.Entropy, na.rm = TRUE), 3),
    Chao1_mean = round(mean(Chao1, na.rm = TRUE), 1),
    Chao1_sd = round(sd(Chao1, na.rm = TRUE), 1),
    Chao1_median = round(median(Chao1, na.rm = TRUE), 1)
  )

print(descriptive_by_health)

write.csv(descriptive_by_phenotype, "Supplementary_Table_descriptive_diversity_by_phenotype.csv", row.names = FALSE)
write.csv(descriptive_by_health, "Supplementary_Table_descriptive_diversity_by_health_status.csv", row.names = FALSE)
