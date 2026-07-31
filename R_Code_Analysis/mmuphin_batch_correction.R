#!/usr/bin/env Rscript

# ==============================================================================
# Batch Correction (MMUPHin) and PERMANOVA Evaluation
#
# Applies MMUPHin batch correction to a species-level abundance matrix using
# study identity (BioProject) as the batch variable, and quantifies the
# proportion of microbial community variance explained by study identity and
# health status (PERMANOVA, Bray-Curtis dissimilarity) before and after
# correction.
#
# Usage:
#   Rscript batch_correct_and_adonis.R <otu_file.csv> <metadata_file> \
#           <corrected_output.tsv> <adonis_output_before.tsv> <adonis_output_after.tsv>
#
# Arguments:
#   otu_file.csv              Species x sample abundance matrix (comma-separated)
#   metadata_file              Sample metadata, including BioProject and
#                               Health_status columns (.tsv = tab-separated,
#                               .csv = comma-separated; detected automatically)
#   corrected_output.tsv       Output path for the batch-corrected abundance matrix
#   adonis_output_before.tsv   Output path for PERMANOVA results before correction
#   adonis_output_after.tsv    Output path for PERMANOVA results after correction
#

library(MMUPHin)
library(vegan)
library(tools)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript batch_correct_and_adonis.R <otu_file.csv> <metadata_file> <corrected_output.tsv> <adonis_output_before.tsv> <adonis_output_after.tsv>\n")
  quit(status = 1)
}

otu_file             <- args[1]
metadata_file        <- args[2]
corrected_output     <- args[3]
adonis_output_before <- args[4]
adonis_output_after  <- args[5]

cat("Loading OTU table...\n")
otu_table <- read.delim(otu_file, sep = ",", row.names = 1, header = TRUE)
cat("OTU table dimensions:", dim(otu_table)[1], "species x", dim(otu_table)[2], "samples\n")

cat("Loading metadata table...\n")
metadata_extension <- tolower(file_ext(metadata_file))
metadata_sep <- if (metadata_extension == "tsv") "\t" else if (metadata_extension == "csv") "," else stop("Unrecognized metadata file extension. Expected .tsv or .csv.")

metadata_table <- read.delim(metadata_file, sep = metadata_sep, row.names = 1, header = TRUE)
cat("Metadata table dimensions:", dim(metadata_table)[1], "samples x", dim(metadata_table)[2], "columns\n")

metadata_table$BioProject <- as.factor(metadata_table$BioProject)

fit_adjust_batch <- adjust_batch(
  feature_abd = otu_table,
  batch = "BioProject",
  covariates = "Health_status",
  data = metadata_table,
  control = list(verbose = FALSE)
)

otu_abd_adj <- fit_adjust_batch$feature_abd_adj

write.table(otu_abd_adj, file = corrected_output, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
cat("Batch correction completed. Corrected OTU table saved to:", corrected_output, "\n")

D_before <- vegdist(t(otu_table), method = "bray")
D_after  <- vegdist(t(otu_abd_adj), method = "bray")

set.seed(1)
fit_adonis_before <- adonis2(D_before ~ BioProject + Health_status, data = metadata_table)
write.table(fit_adonis_before, file = adonis_output_before, sep = "\t", quote = FALSE)
print(fit_adonis_before)

set.seed(1)
fit_adonis_after <- adonis2(D_after ~ BioProject + Health_status, data = metadata_table)
write.table(fit_adonis_after, file = adonis_output_after, sep = "\t", quote = FALSE)
print(fit_adonis_after)

cat("Done. All outputs written successfully.\n")
