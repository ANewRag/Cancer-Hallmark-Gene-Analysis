# Required packages
library("AnnotationDbi")
library("org.Hs.eg.db")
library("DESeq2")

# ==== Merge separated htseq count files into one table ====
tumor_folders = list.files(pattern = "data", full.names = T)
tumor_type = list.files(pattern = "data")
tumor_type = strsplit(tumor_type, "_")
tumor_type = sapply(tumor_type, "[[", 1)

for (t in 1:length(tumor_folders)){
  
  files = list.files(path = tumor_folders[t], pattern = ".htseq.counts.gz$", full.names = T, recursive = T)
  
  data = read.table(files[1], sep="\t")
  table = as.data.frame(matrix(ncol = length(files), nrow = nrow(data)))
  row.names(table) = data[[1]]
  
  for(a in 1:length(files)){
    
    reads = read.table(files[a], sep="\t")
    
    table[[a]] = reads[[2]]
    
    file_name = strsplit(files[a], "/", fixed = T)
    file_name = sapply(file_name, "[[", 3)
    colnames(table)[a] = file_name
    
  }
  
  # ==== Processing of ENSEMBL ids (ENSEMBL id -> GENE SYMBOL) ====
  ens_id = as.character(rownames(table))
  ens_id = strsplit(ens_id, ".", fixed = T)
  ens_id = sapply(ens_id, "[[", 1)
  rownames(table) = ens_id
  rownames(data) = ens_id
  
  data$"V1" = mapIds(org.Hs.eg.db,
                     keys=as.character(rownames(data)),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
  
  # ==== Filter duplicated variants ====
  table_annot = cbind(data[1], table)
  table_annot = table_annot[complete.cases(table_annot), ]
  table_annot = table_annot[!duplicated(as.character(table_annot[[1]])), ] # The first transcript variant was kept!
  rownames(table_annot) = as.character(table_annot[[1]])
  table_annot[[1]] = NULL
  
  # ==== Remove miRNAs from the table ==== 
  table_annot = table_annot[!grepl("^MIR", rownames(table_annot)), ]
  
  # ==== Processing patient's ids ====
  meta = list.files(path = tumor_folders[t], pattern = "gdc_sample_sheet.tsv", full.names = T, recursive = T)
  meta = read.table(meta[1], header = T, sep = "\t", check.names = F)
  tomach = c("01A", "01B", "01C")
  meta = meta[grep(as.character(meta$`Sample ID`), pattern = paste(tomach, collapse = "|")),]
  rownames(meta) = meta$`File ID`
  
  meta_list = rownames(meta)
  int_ids = intersect(colnames(table_annot), meta_list)
  table_annot = table_annot[, int_ids]
  meta = meta[int_ids,]
  
  for (b in 1:ncol(table_annot)) {
    
    colnames(table_annot)[[b]] = as.character(meta[colnames(table_annot)[[b]], "Case ID"])
    
  }
  
  # ==== Removing duplicated samples ==== 
  duplicated_samples = colnames(table_annot[, duplicated(colnames(table_annot))])
  table_annot = table_annot[, !duplicated(colnames(table_annot))] # The first duplicant was kept!!!!!
  # write.table(duplicated_samples, paste(tumor_type[t], "mRNA_suplicated_samples.txt", sep = "_"))
  
  # ==== Normalizing htseq counts data using DESeq ====
  counttable = table_annot
  counttable1 = estimateSizeFactorsForMatrix(counttable)
  normalized_counts = t(t(counttable)/counttable1)
  normalized_counts = round(normalized_counts, digits=0)
  normalized_counts = as.data.frame(normalized_counts)
  
  # ==== Export data ====
  write.table(normalized_counts, paste(tumor_type[t], "mRNA_expression.txt", sep = "_"), sep= "\t", col.names = NA, quote = F)
  
}

# ==== Scaling ====
setwd("~/Documents/TCGA_PanCancer/Expression_tables")

files = list.files(pattern = "mRNA_expression", recursive = T)
tumor_type = strsplit(files, "_")
tumor_type = sapply(tumor_type, "[[", 1)

mrna = read.table(files[1], header = T, sep = "\t", check.names = F, stringsAsFactors = F, row.names = 1)
genes = as.character(rownames(mrna))

for (c in 1:length(files)){
  mrna_c = read.table(files[c], header = T, sep = "\t", check.names = F, stringsAsFactors = F, row.names = 1)
  colmeans_c = colMeans(mrna_c)
  mrna_c = as.data.frame(mapply(function(x,y) round(x/y*1000, digits = 0), x = mrna_c, y = colmeans_c))
  rownames(mrna_c) = genes
  
  write.table(mrna_c, paste(tumor_type[c], "mRNA_scaled_expression.txt", sep = "_"), sep= "\t", col.names = NA, quote = F)
}

