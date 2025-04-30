# 1. Read the sample sheet mapping file
# Adjust the file name/path as necessary
sample_sheet <- read.delim("gdc_sample_sheet.2025-04-18 (3).tsv", 
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# In this sample sheet:
# - The full filename is stored in the 2nd column: sample_sheet[[2]]
# - The corresponding sample ID is stored in the 7th column: sample_sheet[[7]]

# 2. List all RNA-seq count files matching the pattern
files <- list.files(pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$")
if (length(files) == 0) {
  stop("No files found matching the pattern.")
}

# 3. Process the first file to initialize the combined table
first_file <- files[1]
df_first <- read.delim(first_file, header = FALSE, sep = "\t")
df_first <- df_first[-(1:6), ]  # Remove the first 6 rows (metadata/header)
gene_names <- df_first[, 2]      # Gene names assumed in column 2

# Use the entire filename as identifier
file_identifier <- basename(first_file)
sample_id <- sample_sheet[[7]][match(file_identifier, sample_sheet[[2]])]
if (is.na(sample_id)) {
  warning(paste("No matching sample ID found for", file_identifier, "- using the full filename"))
  sample_id <- file_identifier
}

counts_first <- df_first[, 4, drop = FALSE]  # Counts assumed in column 4
colnames(counts_first) <- sample_id

# Create the combined table with gene names and the first sample's counts
combined_table <- data.frame(gene_name = gene_names, counts_first, stringsAsFactors = FALSE, check.names = FALSE)

# 4. Loop over the remaining files and combine count columns
for (file in files[-1]) {
  df <- read.delim(file, header = FALSE, sep = "\t")
  df <- df[-(1:6), ]  # Remove the first 6 rows
  file_identifier <- basename(file)
  sample_id <- sample_sheet[[7]][match(file_identifier, sample_sheet[[2]])]
  if (is.na(sample_id)) {
    warning(paste("No matching sample ID found for", file_identifier, "- using the full filename"))
    sample_id <- file_identifier
  }
  counts <- df[, 4, drop = FALSE]
  colnames(counts) <- sample_id
  combined_table <- cbind(combined_table, counts)
}

# Optionally, write the combined table to a file
write.table(combined_table, file = "combined_gene_counts_with_sample_ids.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 5. Convert count columns to numeric (if necessary)
combined_table[-1] <- lapply(combined_table[-1], function(x) as.numeric(as.character(x)))

# 6. Normalize the counts using DESeq2's size factor estimation
library(DESeq2)

# Create a count matrix (exclude the gene_name column) and set rownames
count_matrix <- as.matrix(combined_table[,-1])
rownames(count_matrix) <- combined_table$gene_name

# Estimate size factors for the count matrix
size_factors <- estimateSizeFactorsForMatrix(count_matrix)

# Normalize the counts by dividing each column by its corresponding size factor
normalized_counts <- t(t(count_matrix) / size_factors)

# Convert the normalized matrix back to a data frame and add gene names as a column
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts <- cbind(gene_name = rownames(normalized_counts), normalized_counts)

# Write the normalized counts table to a file
write.table(normalized_counts, file = "mRNA_scaled_expression.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Optionally, display the first few rows of the normalized counts
head(normalized_counts)
