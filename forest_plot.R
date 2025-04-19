# Required packages
library(survival)
library(survminer)   # for ggsurvplot()
library(forestplot)
library(ggplot2)     # for ggsave()

# ==== SURVIVAL FUNCTION ====
bestcutoff <- function(datavector, clintable) {
  dv_clean <- datavector[!is.na(datavector)]
  breaks <- quantile(dv_clean, probs = seq(0.25, 0.75, by = 0.01), na.rm = TRUE)
  cutoff.table <- t(sapply(breaks, function(z)
    cutoff(datavector = datavector, cutpoint = z, clintable = clintable)))
  colnames(cutoff.table) <- c("cutoff", "pvalue")
  cutoff.table[order(cutoff.table[, 2]), "cutoff"][1]
}
cutoff <- function(datavector, cutpoint, clintable) {
  term <- cut(x = datavector,
              breaks = c(min(datavector, na.rm = TRUE),
                         cutpoint,
                         max(datavector, na.rm = TRUE)),
              labels = FALSE,
              include.lowest = TRUE)
  cox <- summary(coxph(Surv(surv_time, surv_events) ~ term, data = clintable))
  c(cutpoint, cox$sctest[3])
}

# ==== File lists ====
exp_files <- list.files(path = getwd(),
                        pattern = "mRNA_scaled_expression.txt",
                        full.names = TRUE, recursive = TRUE)
clin_file         <- list.files(path = getwd(),
                                pattern = "PanCancer_clinical_table_all190708.txt",
                                full.names = TRUE, recursive = TRUE)
cancer_genes_file <- list.files(path = getwd(),
                                pattern = "Candidate_genes.txt",
                                full.names = TRUE, recursive = TRUE)

# ==== Read clinical & gene‐set tables ====
clinical     <- read.table(clin_file,         sep = "\t", header = TRUE,
                           row.names = 1,    check.names = FALSE)
cancer_genes <- read.table(cancer_genes_file, sep = "\t", header = TRUE,
                           check.names = FALSE)

# ==== Define the 5 hallmarks ====
gs_codes   <- c(
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_E2F_TARGETS"
)
har_labels <- c(
  "DNA_Repair",
  "p53_Pathway",
  "G2M_Checkpoint",
  "Apoptosis",
  "E2F_Targets"
)

# Optional: manual tumor labels (same length/order as exp_files)
manual_tumor_labels <- "Thymoma"
# e.g. manual_tumor_labels <- c("Breast", "Lung", "Colon", ...)

# ==== Loop over all 5 hallmarks ====
for (h in seq_along(gs_codes)) {
  # pick this hallmark
  code          <- gs_codes[h]
  hallmark_label <- har_labels[h]
  hallmark_genes <- unique(
    cancer_genes$gene_symbol[cancer_genes$gs_name == code]
  )
  
  # prepare result tables
  n_tumors      <- length(exp_files)
  forest_input  <- data.frame(matrix(nrow = n_tumors, ncol = 5))
  colnames(forest_input) <- c("Tumor", "p.value", "HR", "CI_Low", "CI_High")
  multi_cox     <- data.frame(matrix(nrow = n_tumors, ncol = 3))
  colnames(multi_cox) <- c("Tumor", "p.value", "HR")
  
  # per‐tumor analysis
  for (i in seq_along(exp_files)) {
    # read expr and decide tumor name
    expr_path <- exp_files[i]
    expr      <- read.table(expr_path, header = TRUE, sep = "\t",
                            check.names = FALSE, row.names = 1)
    if (!is.null(manual_tumor_labels) &&
        length(manual_tumor_labels) == n_tumors) {
      tumor <- manual_tumor_labels[i]
    } else {
      tumor <- strsplit(expr_path, split = "[/,_]")[[1]][9]
    }
    
    # match clinical samples
    samp  <- intersect(colnames(expr), rownames(clinical))
    expr  <- expr[, samp, drop = FALSE]
    clint <- clinical[samp, ]
    
    # survival columns
    surv_time    <- as.numeric(clint[["days_difference"]])
    surv_events  <- as.numeric(clint[["demographic.vital_status"]])
    clint$surv_time   <- surv_time
    clint$surv_events <- surv_events
    
    # hallmark signature
    sig_mat      <- expr[rownames(expr) %in% hallmark_genes, , drop = FALSE]
    hm_signature <- colMeans(sig_mat, na.rm = TRUE)
    
    # best cutoff & grouping
    cutoff_pt    <- as.numeric(bestcutoff(hm_signature, clint))
    exp_category <- factor(
      ifelse(hm_signature >= cutoff_pt, "High", "Low"),
      levels = c("Low", "High")
    )
    
    # univariate Cox
    uni_sum <- summary(coxph(
      Surv(surv_time, surv_events) ~ exp_category,
      data = clint
    ))
    forest_input[i, ] <- list(
      Tumor   = tumor,
      p.value = uni_sum$sctest["pvalue"],
      HR      = round( uni_sum$conf.int[1],   2),
      CI_Low  = round( uni_sum$conf.int[3],   2),
      CI_High = round( uni_sum$conf.int[4],   2)
    )
    
    # multivariate Cox (age, gender, race if present)
    desired_covs <- c(
      "demographic.age_at_index",
      "demographic.gender",
      "demographic.race"
    )
    covs_present <- intersect(make.names(desired_covs), colnames(clint))
    covs_final   <- covs_present[
      colSums(!is.na(clint[, covs_present, drop = FALSE])) > 0
    ]
    rhs  <- paste(c("exp_category", covs_final), collapse = " + ")
    fmla <- as.formula(paste("Surv(surv_time, surv_events) ~", rhs))
    multi_sum <- summary(coxph(fmla, data = clint))
    
    multi_cox[i, ] <- list(
      Tumor   = tumor,
      p.value = multi_sum$sctest["pvalue"],
      HR      = round(multi_sum$conf.int[1,1], 2)
    )
    
    # Kaplan–Meier plot
    fit    <- survfit(Surv(surv_time, surv_events) ~ exp_category,
                      data = clint)
    survplt <- ggsurvplot(
      fit, data = clint, pval = TRUE, risk.table = TRUE,
      conf.int = FALSE, legend.title = "Expression",
      legend.labs = c("Low", "High"),
      palette = c("black", "red"),
      xlab = "Time (months)", ylab = "Survival probability",
      title = paste0(tumor, " - ", hallmark_label),
      risk.table.height = 0.25
    )
    ggsave(
      filename = paste0(tumor, "_", hallmark_label, ".pdf"),
      plot     = survplt$plot,
      width    = 6, height = 5
    )
  }
  
  # write out per‐hallmark result tables
  write.table(
    forest_input,
    file      = paste0(hallmark_label, "_univariate_results.txt"),
    sep       = "\t",
    row.names = FALSE, quote = FALSE
  )
  write.table(
    multi_cox,
    file      = paste0(hallmark_label, "_multivariate_results.txt"),
    sep       = "\t",
    row.names = FALSE, quote = FALSE
  )
}
