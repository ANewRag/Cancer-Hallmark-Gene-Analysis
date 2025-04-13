# ==== LOAD REQUIRED PACKAGES ====
library(survival)
library(ggplot2)

# ==== SURVIVAL FUNCTIONS ====
bestcutoff <- function(datavector, clintable) {
  breaks <- quantile(datavector, probs = seq(0.25, 0.75, by = 0.01))
  cutoff.table <- t(sapply(breaks, function(z) cutoff(datavector = datavector, cutpoint = z, clintable = clintable)))
  colnames(cutoff.table) <- c("cutoff", "pvalue")
  cutoff.table[order(cutoff.table[, 2]), "cutoff"][1]
}

cutoff <- function(datavector, cutpoint, clintable) {
  term <- cut(x = datavector, 
              breaks = c(min(datavector), cutpoint, max(datavector)), 
              labels = FALSE, 
              include.lowest = TRUE)
  cox <- summary(coxph(Surv(surv_time, surv_events) ~ term, data = clintable))
  c(cutpoint, cox$sctest[3])
}

# ==== FILES AND DATA READING ====
exp_files <- list.files(path = getwd(), pattern = "mRNA_scaled_expression.txt", full.names = TRUE, recursive = TRUE)
clin_file <- list.files(path = getwd(), pattern = "PanCancer_clinical_table_all190708.txt", full.names = TRUE, recursive = TRUE)
cancer_genes_file <- list.files(path = getwd(), pattern = "Candidate_genes", full.names = TRUE, recursive = TRUE)

# Read candidate cancer genes (707 unique genes assumed)
cancer_genes <- read.table(cancer_genes_file, sep = "\t", header = TRUE, check.names = FALSE)

# Read clinical table
clinical <- read.table(clin_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

# Select the tumor type (assuming only one is used; here the first file is selected)
a <- 1
expression <- read.table(exp_files[a], header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Extract tumor type from file name (modify index as needed)
tumor_type <- strsplit(exp_files[a], split = "[/,_]")[[1]][9]

# Intersect samples between expression and clinical data and subset accordingly
cases <- intersect(colnames(expression), rownames(clinical))
expression <- expression[, cases]
clintable <- clinical[cases, ]

# Define survival time and events
surv_time <- as.numeric(clintable[["days_difference"]])
surv_events <- as.numeric(clintable[["demographic.vital_status"]])

# ==== PREALLOCATE THE RESULTS MATRIX ====
# We plan to output the following for the gene and 3 covariates (Gender, Race, Age):
# For each variable: p-value, HR, CI lower bound, CI upper bound.
# Adding one column for Gene name gives 1 + (4 variables * 4 values) = 17 columns.
num_genes <- nrow(cancer_genes)
multivariate_surv_analysis <- as.data.frame(matrix(NA, nrow = num_genes, ncol = 17))

colnames(multivariate_surv_analysis) <- c("Gene",
                                          paste(tumor_type, "Gene_P", sep = "_"),
                                          paste(tumor_type, "Gene_HR", sep = "_"),
                                          paste(tumor_type, "Gene_CI_low", sep = "_"),
                                          paste(tumor_type, "Gene_CI_high", sep = "_"),
                                          paste(tumor_type, "Gender_P", sep = "_"),
                                          paste(tumor_type, "Gender_HR", sep = "_"),
                                          paste(tumor_type, "Gender_CI_low", sep = "_"),
                                          paste(tumor_type, "Gender_CI_high", sep = "_"),
                                          paste(tumor_type, "Race_P", sep = "_"),
                                          paste(tumor_type, "Race_HR", sep = "_"),
                                          paste(tumor_type, "Race_CI_low", sep = "_"),
                                          paste(tumor_type, "Race_CI_high", sep = "_"),
                                          paste(tumor_type, "Age_P", sep = "_"),
                                          paste(tumor_type, "Age_HR", sep = "_"),
                                          paste(tumor_type, "Age_CI_low", sep = "_"),
                                          paste(tumor_type, "Age_CI_high", sep = "_")
)

# ==== MAIN LOOP OVER CANDIDATE GENES ====
for (b in 1:nrow(cancer_genes)) {
  tryCatch(
    expr = {
      # Selecting gene expression
      gene_name <- as.character(cancer_genes[b, 1])
      selected_exp <- as.numeric(expression[gene_name, ])
      
      # Determine the best cutoff for dividing patients into high and low expression groups
      cutoff.point <- as.numeric(bestcutoff(datavector = selected_exp, clintable = clintable))
      
      # Divide patients into low/high groups
      exp_category <- ifelse(selected_exp >= cutoff.point, "High", "Low")
      exp_category <- factor(exp_category, levels = c("Low", "High"))
      
      # Create a unified analysis data frame with survival data, exp_category, and covariates.
      analysis_data <- data.frame(
        surv_time = surv_time,
        surv_events = surv_events,
        exp_category = exp_category,
        demographic.gender = clintable[["demographic.gender"]],
        demographic.race = clintable[["demographic.race"]],
        demographic.age_at_index = clintable[["demographic.age_at_index"]]
      )
      
      # Fit the multivariate Cox regression model using the merged data frame
      multivariate_cox_result <- summary(
        coxph(Surv(surv_time, surv_events) ~ exp_category + demographic.gender +
                demographic.race + demographic.age_at_index, data = analysis_data)
      )
      
      # Extract and save results for the gene effect (exp_categoryHigh)
      if ("exp_categoryHigh" %in% rownames(multivariate_cox_result$coefficients)) {
        multivariate_surv_analysis[b, 2] <- round(multivariate_cox_result$coefficients["exp_categoryHigh", "Pr(>|z|)"], 4)
        multivariate_surv_analysis[b, 3] <- round(multivariate_cox_result$conf.int["exp_categoryHigh", "exp(coef)"], 2)
        multivariate_surv_analysis[b, 4] <- round(multivariate_cox_result$conf.int["exp_categoryHigh", "lower .95"], 2)
        multivariate_surv_analysis[b, 5] <- round(multivariate_cox_result$conf.int["exp_categoryHigh", "upper .95"], 2)
      }
      
      # Extract and save results for gender
      if ("demographic.gender" %in% rownames(multivariate_cox_result$coefficients)) {
        multivariate_surv_analysis[b, 6] <- round(multivariate_cox_result$coefficients["demographic.gender", "Pr(>|z|)"], 4)
        multivariate_surv_analysis[b, 7] <- round(multivariate_cox_result$conf.int["demographic.gender", "exp(coef)"], 2)
        multivariate_surv_analysis[b, 8] <- round(multivariate_cox_result$conf.int["demographic.gender", "lower .95"], 2)
        multivariate_surv_analysis[b, 9] <- round(multivariate_cox_result$conf.int["demographic.gender", "upper .95"], 2)
      }
      
      # Extract and save results for race
      if ("demographic.race" %in% rownames(multivariate_cox_result$coefficients)) {
        multivariate_surv_analysis[b, 10] <- round(multivariate_cox_result$coefficients["demographic.race", "Pr(>|z|)"], 4)
        multivariate_surv_analysis[b, 11] <- round(multivariate_cox_result$conf.int["demographic.race", "exp(coef)"], 2)
        multivariate_surv_analysis[b, 12] <- round(multivariate_cox_result$conf.int["demographic.race", "lower .95"], 2)
        multivariate_surv_analysis[b, 13] <- round(multivariate_cox_result$conf.int["demographic.race", "upper .95"], 2)
      }
      
      # Extract and save results for age
      if ("demographic.age_at_index" %in% rownames(multivariate_cox_result$coefficients)) {
        multivariate_surv_analysis[b, 14] <- round(multivariate_cox_result$coefficients["demographic.age_at_index", "Pr(>|z|)"], 4)
        multivariate_surv_analysis[b, 15] <- round(multivariate_cox_result$conf.int["demographic.age_at_index", "exp(coef)"], 2)
        multivariate_surv_analysis[b, 16] <- round(multivariate_cox_result$conf.int["demographic.age_at_index", "lower .95"], 2)
        multivariate_surv_analysis[b, 17] <- round(multivariate_cox_result$conf.int["demographic.age_at_index", "upper .95"], 2)
      }
      
      # Save the gene name
      multivariate_surv_analysis[b, 1] <- gene_name
      
    },
    error = function(e) {
      e_text <- toString(unlist(e))
      error_line <- data.frame(GeneIndex = b, Error = e_text, Time = Sys.time(), stringsAsFactors = FALSE)
      # Optionally, write error_line to a log file
    }
  )
}

# Write out the results table to file
output_filename <- paste0("PanCancer_hallmarkgenes_OS_multivariate_survival_", tumor_type, ".txt")
write.table(multivariate_surv_analysis, output_filename, sep = "\t", quote = FALSE, na = "", col.names = NA)
output_filename <- paste0("PanCancer_hallmarkgenes_OS_multivariate_survival_", tumor_type, ".tsv")
write.table(multivariate_surv_analysis, output_filename, sep = "\t", quote = FALSE, na = "", col.names = NA)
