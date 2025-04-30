# ==== EXTRA PLOTS FOR SURVIVAL ANALYSIS ====
# This script reads the PanCancer survival results and generates:
# 1. Volcano plots of HR vs. p-value for each tumor type
# 2. Heatmap of log2(HR) for the top N significant genes across tumor types (or barplot if only one tumor type)
# 3. Barplot of number of significant genes (FDR < 0.05) per tumor type

# ---- LOAD PACKAGES ----
library(ggplot2)
library(dplyr)
library(pheatmap)

# ---- PARAMETERS ----
results_file <- "PanCancer_hallmarkgenes_OS_survival_results.tsv"
output_dir   <- "plots_extra"

top_n_genes  <- 50      # Number of top genes for heatmap/barplot
alpha_sig    <- 0.05    # FDR threshold for significance

# ---- READ DATA ----
res <- read.delim(results_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(res) > 0 && colnames(res)[1] == "") colnames(res)[1] <- "gene"
res$gene <- as.character(res$gene)

# Identify p-value and HR columns by their suffix
p_cols  <- grep("_pvalue$", names(res), value = TRUE)
hr_cols <- grep("_HR$", names(res), value = TRUE)

if (length(p_cols) == 0) stop("No p-value columns detected. They should end with '_pvalue'.")
if (length(hr_cols) == 0) stop("No HR columns detected. They should end with '_HR'.")

# Create output directory
if (!dir.exists(output_dir)) dir.create(output_dir)

# ---- 1. Volcano plots per tumor type ----
for (pcol in p_cols) {
  hr_col <- sub("_pvalue$", "_HR", pcol)
  if (!(hr_col %in% names(res))) next
  d <- data.frame(
    gene = res$gene,
    p    = as.numeric(res[[pcol]]),
    HR   = as.numeric(res[[hr_col]])
  ) %>% mutate(
    log2HR  = log2(HR),
    negLogP = -log10(p)
  )
  p <- ggplot(d, aes(x = log2HR, y = negLogP)) +
    geom_point(alpha = 0.5, na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = paste0("Volcano plot: ", sub("_pvalue$", "", pcol)), x = "log2(HR)", y = "-log10(p-value)") +
    theme_minimal()
  ggsave(file.path(output_dir, paste0(sub("_pvalue$", "", pcol), "_volcano.png")), p, width = 6, height = 5)
}

# ---- 2. Heatmap or barplot of top N significant genes ----
# Compute FDR-adjusted p-values per tumor type
for (pcol in p_cols) {
  res[[paste0(pcol, "_adj")]] <- p.adjust(as.numeric(res[[pcol]]), method = "fdr")
}
adj_cols <- grep("_pvalue_adj$", names(res), value = TRUE)
res$min_adj <- if(length(adj_cols)==1) res[[adj_cols]] else apply(res[, adj_cols, drop=FALSE], 1, min, na.rm=TRUE)

top_genes <- res %>% arrange(min_adj) %>% slice(1:top_n_genes) %>% pull(gene)
heat_df <- res %>% filter(gene %in% top_genes) %>% select(gene, all_of(hr_cols))
mat <- as.matrix(heat_df[, hr_cols, drop=FALSE])
rownames(mat) <- heat_df$gene
mat_log2 <- log2(mat)

if (ncol(mat_log2) >= 2 && nrow(mat_log2) >= 2) {
  # clustered heatmap
  pheatmap(mat_log2, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=TRUE, fontsize_row=6,
           main=paste0("Top ", top_n_genes, " Genes by FDR"), filename=file.path(output_dir, "heatmap_top_genes.png"))
} else if (ncol(mat_log2)==1 && nrow(mat_log2) >= 2) {
  # single tumor type: horizontal barplot of log2(HR)
  df_single <- data.frame(gene=rownames(mat_log2), log2HR=mat_log2[,1])
  tumor <- colnames(mat_log2)
  p1 <- ggplot(df_single, aes(x=reorder(gene, log2HR), y=log2HR)) +
    geom_col() + coord_flip() +
    labs(title=paste0("Log2(HR) for Top ", top_n_genes, " Genes: ", tumor), x="Gene", y="log2(HR)") +
    theme_bw()
  ggsave(file.path(output_dir, paste0(tumor, "_top_genes_barplot.png")), p1, width=6, height=top_n_genes*0.2 + 2)
} else {
  warning("Not enough data for heatmap or barplot: got ", nrow(mat_log2), " genes and ", ncol(mat_log2), " tumor types.")
}

# ---- 3. Barplot of significant gene counts per tumor ----
sig_counts <- data.frame(tumor=sub("_pvalue_adj$", "", adj_cols), count=sapply(adj_cols, function(col) sum(res[[col]] < alpha_sig, na.rm=TRUE)))
p_bar <- ggplot(sig_counts, aes(x=reorder(tumor, -count), y=count)) + geom_col() +
  labs(title="Number of Significant Genes per Tumor Type", x="Tumor Type", y="Count (FDR < 0.05)") +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(output_dir, "sig_gene_counts.png"), p_bar, width=8, height=5)

message("All extra plots saved in ", output_dir)
