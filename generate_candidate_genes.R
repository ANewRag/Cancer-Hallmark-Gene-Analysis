# =============================
# Generate Candidate_genes.txt (3 columns)
# =============================

# Install and load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("msigdbr")) BiocManager::install("msigdbr")
if (!require("dplyr")) install.packages("dplyr")

library(msigdbr)
library(dplyr)

# Settings
species <- "Homo sapiens"
category <- "H"
selected_sets <- c("HALLMARK_APOPTOSIS",
                   "HALLMARK_P53_PATHWAY",
                   "HALLMARK_DNA_REPAIR",
                   "HALLMARK_E2F_TARGETS",
                   "HALLMARK_G2M_CHECKPOINT")

# Get MSigDB genes
hallmark_genes <- msigdbr(species = species, category = category)

# Filter and select 3 columns
candidate_genes <- hallmark_genes %>%
  filter(gs_name %in% selected_sets) %>%
  select(gene_symbol, gs_name, gs_description) %>%
  distinct() %>%
  arrange(gene_symbol)

# Save to file
write.table(candidate_genes,
            file = "Candidate_genes.txt",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

cat("✅ Candidate_genes.txt created with", nrow(candidate_genes), "rows and 3 columns.\n")

