##############################
# LOAD REQUIRED PACKAGES
##############################
library(survival)
library(forestplot)
library(survminer)     # Using survminer for modern survival plots

##############################
# DEBUGGING AND LOGGING SETUP
##############################
debug_mode <- TRUE
log_debug <- function(message) {
  if (debug_mode) {
    cat(message, "\n")
  }
}

##############################
# SURVIVAL FUNCTIONS
##############################
bestcutoff <- function(datavector, clintable) {
  # Calculate candidate cutoff values using quantiles
  breaks <- quantile(datavector, probs = seq(0.25, 0.75, by = 0.01), na.rm = TRUE)
  cutoff.table <- t(sapply(breaks, function(z) cutoff(datavector = datavector, cutpoint = z, clintable = clintable)))
  colnames(cutoff.table) <- c("cutoff", "pvalue")
  
  log_debug(paste("Cutoff table:\n", paste(capture.output(print(cutoff.table)), collapse = "\n")))
  best_index <- which.min(cutoff.table[, 2])
  best_cutoff <- cutoff.table[best_index, "cutoff"]
  log_debug(paste("Selected best cutoff:", best_cutoff, "with minimum p-value:", cutoff.table[best_index, "pvalue"]))
  best_cutoff
}

cutoff <- function(datavector, cutpoint, clintable) {
  # Create a categorical variable using the candidate cutoff value
  term <- cut(x = datavector, 
              breaks = c(min(datavector, na.rm = TRUE), cutpoint, max(datavector, na.rm = TRUE)),
              labels = FALSE, include.lowest = TRUE)
  log_debug(paste("Testing cutpoint:", round(cutpoint, 3), 
                  "Group sizes:", paste(names(table(term)), table(term), collapse = " ")))
  cox <- withCallingHandlers({
    summary(coxph(Surv(surv_time, surv_events) ~ term, data = clintable))
  }, warning = function(w) {
    log_debug(paste("Warning for cutoff", round(cutpoint,3), ":", conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
  log_debug(paste("Cutpoint:", round(cutpoint, 3), "score test p-value:", cox$sctest[3]))
  c(cutpoint, cox$sctest[3])
}

##############################
# INITIALIZE RESULT DATAFRAMES
##############################
exp_files <- list.files(path = getwd(), pattern = "mRNA_scaled_expression.txt", full.names = TRUE, recursive = TRUE)
num_files <- length(exp_files)

forest_input <- as.data.frame(matrix(nrow = num_files, ncol = 5))
colnames(forest_input) <- c("Tumor", "p.value", "HR", "CI_Low", "CI_High")

multi_cox <- as.data.frame(matrix(nrow = num_files, ncol = 3))
colnames(multi_cox) <- c("Tumor", "p.value", "HR")

##############################
# FILES AND DATA READING
##############################
clin_file <- list.files(path = getwd(), pattern = "PanCancer_clinical_table_all190708.txt", full.names = TRUE, recursive = TRUE)
cancer_genes_file <- list.files(path = getwd(), pattern = "Candidate_genes.txt", full.names = TRUE, recursive = TRUE)

cancer_genes <- read.table(cancer_genes_file, sep = "\t", header = TRUE, check.names = FALSE)
clinical <- read.table(clin_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

##############################
# SELECT MAIN CANCER HALLMARK
##############################
selected_hallmark <- 6  # Change this number as desired
if (selected_hallmark == 1){
  hallmark <- "Sustaining proliferative signaling"
  hallmark_genes <- unique(as.character(cancer_genes[which(cancer_genes$Hallmark_feature == "oncogenes"), 1]))
} else if (selected_hallmark == 2){
  hallmark <- "Evading growth suppressors"
  hallmark_genes <- unique(as.character(cancer_genes[which(cancer_genes$Hallmark_feature == "tumor_suppressor_genes"), 1]))
} else if (selected_hallmark == 3){
  hallmark <- "Inducing angiogenesis"
  hallmark_genes <- unique(as.character(cancer_genes[which(cancer_genes$Main_hallmark_feature == "inducing_angiogenesis"), 1]))
} else if (selected_hallmark == 4){
  hallmark <- "Genome instability"
  hallmark_genes <- unique(as.character(cancer_genes[which(cancer_genes$Main_hallmark_feature == "genome_instability"), 1]))
} else if (selected_hallmark == 5){
  hallmark <- "Deregulation of cellular energetics"
  hallmark_genes <- unique(as.character(cancer_genes[which(cancer_genes$Main_hallmark_feature == "deregulation_of_cellular_energetics"), 1]))
} else if (selected_hallmark == 6){
  hallmark <- "Activating invasion metastasis"
  hallmark_genes <- unique(as.character(cancer_genes[which(cancer_genes$Main_hallmark_feature %in% c("activating_invasion_metastasis", 
                                                                                                     "activating_invasion_metastasis/cell_motility", 
                                                                                                     "activating_invasion_metastasis/rho_family_ GTPases", 
                                                                                                     "activating_invasion_metastasis/invasion")), 1]))
} else if (selected_hallmark == 7){
  hallmark <- "Resisting cell death"
  hallmark_genes <- unique(as.character(cancer_genes[which(cancer_genes$Main_hallmark_feature %in% c("resisting_cell_death/apoptotic_pathway_members", 
                                                                                                     "resisting_cell_death/TP53_pathway")), 1]))
} else if (selected_hallmark == 8){
  hallmark <- "Enabling replicative immortality"
  hallmark_genes <- unique(as.character(cancer_genes[which(cancer_genes$Main_hallmark_feature %in% c("enabling_replicative_immortality", 
                                                                                                     "enabling_replicative_immortality/telomerase_activity", 
                                                                                                     "enabling_replicative_immortality/senescence_secretome")), 1]))
}

##############################
# SURVIVAL ANALYSIS PER HALLMARK
##############################
for (a in 1:length(exp_files)) {
  log_debug(paste("Processing expression file:", exp_files[a]))
  expression <- read.table(exp_files[a], header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  
  # Extract tumor type (adjust the index as needed)
  tumor_parts <- strsplit(exp_files[a], split = "[/,_]")[[1]]
  tumor_type <- tumor_parts[9]
  
  cases <- intersect(colnames(expression), rownames(clinical))
  expression <- expression[, cases, drop = FALSE]
  clintable <- clinical[cases, ]
  
  # Survival data (assumed: column 3 = time; column 4 = event)
  surv_time <- as.numeric(clintable[["days_difference"]])
  surv_events <- as.numeric(clintable[["demographic.vital_status"]])
  
  # Compute hallmark signature (mean of available hallmark genes)
  common_genes <- intersect(hallmark_genes, rownames(expression))
  if (length(common_genes) == 0) {
    log_debug(paste("No hallmark genes found in expression data for tumor type:", tumor_type))
    next
  }
  expression_sign <- expression[common_genes, , drop = FALSE]
  hm_signature <- colMeans(expression_sign, na.rm = TRUE)
  
  # Determine best cutoff for the hallmark signature
  cutoff.point <- bestcutoff(datavector = hm_signature, clintable = clintable)
  log_debug(paste("Best cutoff for hallmark signature in tumor", tumor_type, ":", cutoff.point))
  
  # Divide patients into "High" and "Low" groups (vectorized)
  exp_category <- ifelse(hm_signature >= cutoff.point, "High", "Low")
  exp_category <- factor(exp_category, levels = c("Low", "High"))
  
  # Univariate Cox regression
  cox_result <- withCallingHandlers({
    coxph(Surv(surv_time, surv_events) ~ exp_category, data = clintable)
  }, warning = function(w) {
    log_debug(paste("Coxph warning for tumor", tumor_type, ":", conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
  sum_cox <- summary(cox_result)
  
  forest_input[a, "Tumor"]   <- tumor_type
  forest_input[a, "p.value"] <- as.numeric(sum_cox$sctest['pvalue'])
  forest_input[a, "HR"]      <- as.numeric(round(sum_cox$conf.int[1], digits = 2))
  forest_input[a, "CI_Low"]  <- as.numeric(round(sum_cox$conf.int[3], digits = 2))
  forest_input[a, "CI_High"] <- as.numeric(round(sum_cox$conf.int[4], digits = 2))
  
  # Multivariate Cox regression (including selected clinical covariates)
  original_clincols <- colnames(clintable)
  colnames(clintable) <- make.names(colnames(clintable))
  clin_names <- make.names(c("Stage", "Gender(0=female;1=male)", 
                             "Race(1=white;2=asian;3=black/african american)", 
                             "Grade(0(LowGrade)=G1+G2;1(HighGrade)=G3+G4)", "Age"))
  clin_names_allna <- sapply(clintable[, clin_names, drop = FALSE], function(x) !all(is.na(x)))
  clin_names_nona <- names(clin_names_allna[clin_names_allna == TRUE])
  
  if (length(clin_names_nona) > 0) {
    formula_str <- paste("Surv(surv_time, surv_events) ~ exp_category +", paste(clin_names_nona, collapse = " + "))
    fmla <- as.formula(formula_str)
    multi_cox_result <- withCallingHandlers({
      coxph(fmla, data = clintable)
    }, warning = function(w) {
      log_debug(paste("Multivariate Coxph warning for tumor", tumor_type, ":", conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
    sum_multi_cox <- summary(multi_cox_result)
    
    multi_cox[a, "Tumor"]   <- tumor_type
    multi_cox[a, "p.value"] <- as.numeric(sum_multi_cox$sctest['pvalue'])
    multi_cox[a, "HR"]      <- as.numeric(round(sum_multi_cox$conf.int[1, 1], digits = 2))
  } else {
    log_debug(paste("No available clinical covariates for multivariate analysis in tumor:", tumor_type))
  }
  
  # --- Survival Plot using survminer::ggsurvplot ---
  # Uncomment the following block to generate survival plots for the hallmark signature
  #
  # fit <- survfit(Surv(surv_time, surv_events) ~ exp_category, data = clintable)
  # pdf(file = paste0(tumor_type, "_", hallmark, "_SurvivalPlot.pdf"))
  # p <- ggsurvplot(fit, data = clintable,
  #                 pval = TRUE,
  #                 risk.table = TRUE,
  #                 legend.title = "Expression",
  #                 legend.labs = c("Low", "High"),
  #                 xlab = "Time (months)",
  #                 ylab = "Survival probability",
  #                 title = paste(tumor_type, "-", hallmark))
  # print(p$plot)
  # dev.off()
}

# Save multivariate Cox results to file
write.table(multi_cox, file = paste(hallmark, "multivariate_results.txt", sep = "_"),
            sep = "\t", col.names = NA, na = "")

##############################
# SURVIVAL ANALYSIS FOR TUMOR MUTATION BURDEN (TMB)
##############################
tumor_types <- sapply(strsplit(exp_files, split = "[/,_]"), function(x) x[9])
for (a in 1:length(tumor_types)) {
  tumor <- tumor_types[a]
  log_debug(paste("Processing TMB analysis for tumor:", tumor))
  
  clinical_tumor <- clinical[clinical$Project_ID == tumor, ]
  if (nrow(clinical_tumor) == 0) {
    log_debug(paste("No clinical data for tumor:", tumor))
    next
  }
  
  mut_burd <- clinical_tumor[!is.na(clinical_tumor$Total_Number_Of_Mutation), ]
  if (nrow(mut_burd) == 0) {
    log_debug(paste("No mutation burden data for tumor:", tumor))
    next
  }
  
  mut_burd_values <- as.numeric(mut_burd$Total_Number_Of_Mutation)
  surv_time <- as.numeric(mut_burd[[3]])
  surv_events <- as.numeric(mut_burd[[4]])
  
  cutoff.point <- bestcutoff(datavector = mut_burd_values, clintable = mut_burd)
  log_debug(paste("Best TMB cutoff for tumor", tumor, ":", cutoff.point))
  
  exp_category <- ifelse(mut_burd_values >= cutoff.point, "High", "Low")
  exp_category <- factor(exp_category, levels = c("Low", "High"))
  
  cox_result <- withCallingHandlers({
    coxph(Surv(surv_time, surv_events) ~ exp_category, data = mut_burd)
  }, warning = function(w) {
    log_debug(paste("TMB Coxph warning for tumor", tumor, ":", conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
  sum_cox <- summary(cox_result)
  
  idx <- which(forest_input$Tumor == tumor)
  if (length(idx) == 0) { idx <- a }
  forest_input[idx, "Tumor"]   <- tumor
  forest_input[idx, "p.value"] <- as.numeric(sum_cox$sctest['pvalue'])
  forest_input[idx, "HR"]      <- as.numeric(round(sum_cox$conf.int[1], digits = 2))
  forest_input[idx, "CI_Low"]  <- as.numeric(round(sum_cox$conf.int[3], digits = 2))
  forest_input[idx, "CI_High"] <- as.numeric(round(sum_cox$conf.int[4], digits = 2))
  
  # --- TMB Survival Plot using survminer::ggsurvplot ---
  # Uncomment to generate TMB survival plots
  #
  # fit_tmb <- survfit(Surv(surv_time, surv_events) ~ exp_category, data = mut_burd)
  # pdf(file = paste0(tumor, "_TMB_SurvivalPlot.pdf"))
  # p_tmb <- ggsurvplot(fit_tmb, data = mut_burd,
  #                     pval = TRUE,
  #                     risk.table = TRUE,
  #                     legend.title = "TMB Group",
  #                     legend.labs = c("Low", "High"),
  #                     xlab = "Time (months)",
  #                     ylab = "Survival probability",
  #                     title = paste("TMB Survival for", tumor))
  # print(p_tmb$plot)
  # dev.off()
}

# Sort and format forest_input for plotting
forest_input_sorted <- forest_input[order(as.numeric(forest_input$p.value)), ]
forest_input_sorted$p.value <- format(as.numeric(forest_input_sorted$p.value), digits = 2, scientific = TRUE)

##############################
# CREATE FOREST PLOT FOR TMB RESULTS
##############################
header <- as.data.frame(matrix(nrow = 1, ncol = 3))
header[1, ] <- c("Tumor", "HR (CI)", "p-values")
colnames(header) <- c("Tumor", "HR (CI)", "p-values")

text <- as.data.frame(matrix(nrow = nrow(forest_input_sorted), ncol = 3))
colnames(text) <- c("Tumor", "HR (CI)", "p-values")
text[, "Tumor"]   <- forest_input_sorted$Tumor
text[, "HR (CI)"] <- paste(forest_input_sorted$HR, paste("(", forest_input_sorted$CI_Low, "-", forest_input_sorted$CI_High, ")", sep = ""), sep = " ")
text[, "p-values"] <- forest_input_sorted$p.value

plot_supp_info <- rbind(as.matrix(header), as.matrix(text))

plot_data <- forest_input_sorted[, c("HR", "CI_Low", "CI_High")]
colnames(plot_data) <- c("HR", "CI-low", "CI-high")
header2 <- as.data.frame(matrix(nrow = 1, ncol = 3))
colnames(header2) <- c("HR", "CI-low", "CI-high")
plot_input_data <- rbind(header2, plot_data)

pdf("TMB_plot.pdf")
forestplot(plot_supp_info,
           plot_input_data,
           title = "Tumor Mutation Burden",
           hrzl_lines = list("2" = gpar(lty = 1)),
           clip = c(0, 4.5),
           new_page = FALSE,
           zero = 1,
           lwd.zero = 0.8,
           xlog = FALSE,
           graph.pos = 2,
           xlab = "HR-values",
           colgap = unit(1, "cm"),
           col = fpColors(box = "black", lines = "black", zero = "gray50"),
           cex = 0.9,
           lineheight = "auto",
           boxsize = 0.25,
           lwd.ci = 2,
           ci.vertices = TRUE,
           ci.vertices.height = 0.1,
           fn.ci_norm = "fpDrawNormalCI",
           txt_gp = fpTxtGp(cex = 0.8, ticks = gpar(fontfamily = "", cex = 0.9),
                            xlab = gpar(fontfamily = "", cex = 1))
)
dev.off()
