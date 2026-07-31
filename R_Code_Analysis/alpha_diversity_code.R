#!/usr/bin/env Rscript

# ==============================================================================
# Alpha Diversity Estimation 
# Run in job submission environment (HTCondor / ReCaS-Bari cluster)
# Output further processed and visualized in RStudio (see figures/ scripts)
# ==============================================================================

library(vegan)
library(iNEXT)
library(Tjazi)

input_file <- "/lustrehome/babluuniba2022/bablu/condor/Data_Analysis/JUNE_2025/otu_table_7452.tsv"
output_file <- "/lustrehome/babluuniba2022/bablu/condor/Data_Analysis/JUNE_2025/alpha_diversity_7452_species.tsv"

otus <- read.delim(input_file, sep = "\t", row.names = 1, header = TRUE)

cat("Dimensions of OTU table:", dim(otus), "\n")

otus_table <- apply(otus, c(1, 2), function(x) as.numeric(as.character(x)))

print(head(otus_table))

# --- Alpha Diversity Calculation ---
alpha_diversity <- get_asymptotic_alpha(species = otus_table, verbose = FALSE)
print(head(alpha_diversity))

write.table(alpha_diversity, output_file, sep = "\t", row.names = TRUE, quote = FALSE)

cat("Alpha diversity results saved to:", output_file, "\n")
