setwd("/Users/bablu/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/CODE/Thesis_Figures_R")

library(ggplot2)
library(dplyr)

performance <- read.delim("per_disease_performance_healthy_vs_disease_SVM_RBF_PFI_20pct.tsv")

performance <- data.frame(
  Disease   = c("CRC","PD","Obesity","CD","UC","AP","ACVD","LV","AS","GC","AMD","GDM"),
  N_disease = c(205, 157, 113, 90, 58, 42, 34, 29, 21, 20, 19, 10),
  AUC       = c(0.979362, 0.952096, 0.918362, 0.971829, 0.963776, 0.954820,
                0.939903, 0.953774, 0.955473, 0.958297, 0.857295, 0.948052)
)

corr_test <- cor.test(performance$N_disease, performance$AUC, method = "spearman")
corr_auc  <- round(corr_test$estimate, 2)
pval_auc  <- round(corr_test$p.value, 2)
sig_label <- ifelse(corr_test$p.value < 0.05, "significant", "not significant")

median_auc <- median(performance$AUC)

p <- ggplot(performance, aes(x = N_disease, y = AUC)) +
  geom_smooth(
    method = "lm", se = TRUE, level = 0.95,
    color = "#B22222", linetype = "dashed", linewidth = 0.9, alpha = 0.15
  ) +
  geom_point(
    size = 5, shape = 21, fill = "#2E6F95", color = "black", stroke = 1.1, alpha = 0.9
  ) +
  geom_text(
    aes(label = paste0(Disease, " (n=", N_disease, ")")),
    hjust = -0.12, vjust = -0.6, size = 3.4, fontface = "plain", color = "#1a1a1a"
  ) +
  geom_hline(yintercept = median_auc, linetype = "dotted", color = "gray40", linewidth = 0.5, alpha = 0.6) +
  annotate(
    "text",
    x = max(performance$N_disease) * 0.98, y = median_auc + 0.003,
    label = paste0("Median AUC = ", round(median_auc, 3)),
    hjust = 1, vjust = 0, size = 3.2, fontface = "italic", color = "gray40"
  ) +
  scale_y_continuous(limits = c(min(performance$AUC) - 0.03, 1.0)) +
  labs(
    x = "Disease Sample Size (N)",
    y = "AUC (Healthy vs Disease)",
    title = "Model Discriminative Power",
    subtitle = paste0("SVM_RBF_PFI — Spearman r = ", corr_auc, ", p = ", pval_auc, " (", sig_label, ")")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(linetype = "dashed", linewidth = 0.3, color = "grey85"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 11, hjust = 0),
    plot.margin = margin(10, 25, 10, 10)
  )

print(p)

ggsave("AUC_vs_disease_sample_size.png", p, width = 9, height = 6.5, dpi = 600)
ggsave("AUC_vs_disease_sample_size.pdf", p, width = 9, height = 6.5)
