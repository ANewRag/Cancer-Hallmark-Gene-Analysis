# ==== LOAD REQUIRED PACKAGES ====
library(survival)
library(beeswarm)
library(ggplot2)

# ==== DEBUGGING AND LOGGING FUNCTIONS ====
# Optionally set a debug flag to enable/disable verbose logging:
debug_mode <- TRUE

log_debug <- function(message) {
  if(debug_mode) {
    cat(message, "\n")
  }
}

# ==== SURVIVAL FUNCTIONS ====
bestcutoff <- function(datavector, clintable) {
  # Calculate a grid of candidate cutoff values
  breaks <- quantile(datavector, probs = seq(0.25, 0.75, by = 0.01), na.rm = TRUE)
  
  # Evaluate each cutoff candidate with the cutoff() function
  cutoff.table <- t(sapply(breaks, function(z) cutoff(datavector = datavector, cutpoint = z, clintable = clintable)))
  colnames(cutoff.table) <- c("cutoff", "pvalue")
  
  log_debug(paste("Cutoff table for current gene:\n", paste(capture.output(print(cutoff.table)), collapse="\n")))
  
  best.index <- which.min(cutoff.table[,2])
  best.cut <- cutoff.table[best.index, "cutoff"]
  log_debug(paste("Selected best cutoff:", best.cut, "with minimum p-value:", cutoff.table[best.index, "pvalue"]))
  best.cut
}

cutoff <- function(datavector, cutpoint, clintable) {
  # Create a categorical variable using the candidate cutoff
  term <- cut(x = datavector, breaks = c(min(datavector), cutpoint, max(datavector)),
              labels = FALSE, include.lowest = TRUE)
  
  # Log the group sizes for this cutoff candidate
  log_debug(paste("Testing cutpoint:", round(cutpoint, 3), "Group sizes:", paste(names(table(term)), table(term), collapse=" ")))
  
  # Fit the cox model with warning handling
  cox <- withCallingHandlers({
    summary(coxph(Surv(surv_time, surv_events) ~ term, data = clintable))
  },
  warning = function(w) {
    log_debug(paste("Warning for cutoff", round(cutpoint,3), ":", conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
  
  log_debug(paste("Cutpoint:", round(cutpoint, 3), "score test p-value:", cox$sctest[3]))
  c(cutpoint, cox$sctest[3])
}

# ==== FILES AND DATA READING ====
# List files for expression, clinical, and candidate gene lists
exp_files <- list.files(path = getwd(), pattern = "mRNA_scaled_expression.txt", full.names = TRUE, recursive = TRUE)
clin_file <- list.files(path = getwd(), pattern = "PanCancer_clinical_table_all190708.txt", full.names = TRUE, recursive = TRUE)
cancer_genes_file <- list.files(path = getwd(), pattern = "Candidate_genes.txt", full.names = TRUE, recursive = TRUE)

# Read candidate cancer genes (707 unique genes)
cancer_genes <- read.table(cancer_genes_file, sep = "\t", header = TRUE, check.names = FALSE)

# Read clinical table
clinical <- read.table(clin_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

# Initialize a results data frame (adjust dimensions accordingly)
surv_analysis <- as.data.frame(matrix(nrow = 911, ncol = 107))

# ==== MAIN LOOP OVER EXPRESSION FILES ====
for (a in 1:length(exp_files)) {
  
  # Read the expression table
  expression <- read.table(exp_files[a], header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  log_debug(paste("Processing file:", exp_files[a]))
  log_debug(paste("All values NA?", all(is.na(expression))))
  
  tumor_type <- strsplit(exp_files[a], split = "[/,_]")[[1]]
  # You may need to adapt this line to extract the desired tumor type string from tumor_type vector
  
  # Intersect clinical and expression samples
  cases <- intersect(colnames(expression), rownames(clinical))
  expression <- expression[, cases]
  log_debug(paste("Post-intersection all values NA?", all(is.na(expression))))
  
  # Subset clinical table for the current cases
  clintable <- clinical[cases, ]
  
  # Define survival time and events (ensure these column names are correct)
  surv_time <- as.numeric(clintable[["days_difference"]])
  surv_events <- as.numeric(clintable[["demographic.vital_status"]])
  
  # Loop over candidate genes
  for (b in 1:nrow(cancer_genes)) {
    tryCatch({
      # Optionally, enable interactive debugging for a specific gene:
      # if(cancer_genes[b,1] == "BAX") browser()
      
      # ==== Select gene expression ====
      gene_name <- as.character(cancer_genes[b, 1])
      selected_exp <- as.numeric(expression[gene_name, ])
      log_debug(paste("Processing gene:", gene_name))
      log_debug("Expression summary:")
      print(summary(selected_exp))
      
      # Optionally, save a histogram (uncomment if needed)
      # pdf(paste0("hist_", gene_name, ".pdf"))
      # hist(selected_exp, main = paste("Histogram for", gene_name), xlab = "Expression")
      # dev.off()
      
      # ==== Determine the best cutoff ====
      cutoff.point <- as.numeric(bestcutoff(datavector = selected_exp, clintable = clintable))
      log_debug(paste("Best cutoff for gene", gene_name, ":", cutoff.point))
      
      # ==== Divide patients into 'Low' and 'High' expression groups ====
      exp_category <- ifelse(selected_exp >= cutoff.point, "High", "Low")
      exp_category <- factor(exp_category, levels = c("Low", "High"))
      log_debug("Group sizes after categorization:")
      print(table(exp_category))
      log_debug(paste("Events in Low group:", sum(surv_events[exp_category == "Low"], na.rm = TRUE)))
      log_debug(paste("Events in High group:", sum(surv_events[exp_category == "High"], na.rm = TRUE)))
      
      # ==== Fit Cox Regression Model with increased iteration limit ====
      cox_result <- withCallingHandlers({
        coxph(Surv(surv_time, surv_events) ~ exp_category,
              control = coxph.control(iter.max = 100))
      }, warning = function(w) {
        log_debug(paste("Coxph warning for gene", gene_name, ":", conditionMessage(w)))
        invokeRestart("muffleWarning")
      })
      log_debug(paste("Cox model summary for gene:", gene_name))
      print(summary(cox_result))
      
      # ==== Save model results ====
      surv_analysis[b, a*4]     <- as.numeric(summary(cox_result)$sctest['pvalue'])
      colnames(surv_analysis)[a*4] <- paste(tumor_type[1], "pvalue", sep = "_")
      
      surv_analysis[b, a*4 + 1] <- as.numeric(round(summary(cox_result)$conf.int[1], digits = 2))
      colnames(surv_analysis)[a*4 + 1] <- paste(tumor_type[1], "HR", sep = "_")
      
      surv_analysis[b, a*4 + 2] <- as.numeric(round(summary(cox_result)$conf.int[3], digits = 2))
      colnames(surv_analysis)[a*4 + 2] <- paste(tumor_type[1], "CI_low", sep = "_")
      
      surv_analysis[b, a*4 + 3] <- as.numeric(round(summary(cox_result)$conf.int[4], digits = 2))
      colnames(surv_analysis)[a*4 + 3] <- paste(tumor_type[1], "CI_high", sep = "_")
      
      # ==== Create and print Beeswarm plot ====
      exp_category_bee <- cut(selected_exp, breaks = quantile(selected_exp, c(0, 0.25, 0.75, 1)),
                              labels = c("low", "mid", "low"), include.lowest = TRUE)
      plot_obj <- ggplot(mapping = aes(x = "BAX", y = selected_exp, colour = exp_category_bee)) +
        geom_jitter(width = 0.5, show.legend = FALSE, size = 3) +
        scale_color_manual(values = c("#999999", "black")) +
        labs(y = "Expression", title = gene_name) +
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.line = element_line(colour = "black"),
              plot.title = element_text(size = 20, hjust = 0.5))
      print(plot_obj)
      
    }, error = function(e) {
      log_debug(paste("Error for gene", as.character(cancer_genes[b, 1]), "at iteration", b, ":", conditionMessage(e)))
      error_line <- data.frame(gene = cancer_genes[b, 1],
                               error_message = toString(e),
                               time = Sys.time(),
                               stringsAsFactors = FALSE)
      # Optionally, write the error log to a file:
      # write.table(error_line, file = paste0("Error_Warnings/", paste(tumor_type[1], "error_pancancer.txt", sep = "_")),
      #             append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    })
  }
  
  log_debug(paste("Finished processing tumor type: ", tumor_type[1]))
}

# ==== SAVE RESULTS ====
# Add gene names (if needed) before writing
surv_analysis[[1]] <- as.character(cancer_genes[[1]])
surv_analysis[[2]] <- as.character(cancer_genes[[2]])
surv_analysis[[3]] <- as.character(cancer_genes[[3]])
write.table(surv_analysis, "PanCancer_hallmarkgenes_OS_survival_results.txt",
            sep = "\t", quote = FALSE, na = "", col.names = NA)

write.table(surv_analysis, "PanCancer_hallmarkgenes_OS_survival_results.tsv",
            sep = "\t", quote = FALSE, na = "", col.names = NA)

#############################
# MULTIPLE TESTING CORRECTION FOR HALLMARK GENES
#############################
#Need to maybe add a task for bonferroni correction
