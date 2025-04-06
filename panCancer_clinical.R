library(data.table)
sample_sheet <- read.delim("gdc_sample_sheet.2025-03-27.tsv", 
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)

clinical_sheet <- fread("clinical.tsv")

removed_dup <- clinical_sheet[!duplicated(clinical_sheet[[2]]), ]

selected_columns <- removed_dup[, .SD, .SDcols = c(10, 17, 18, 29)]

setDT(sample_sheet)

# Add a new column 'sample_id' by matching the case id from selected_columns to sample_sheet$case_id
selected_columns[, sample_id := sample_sheet[match(get(names(selected_columns)[1]), Case.ID), Sample.ID]]

# Reorder the columns to insert 'sample_id' as the second column
current_names <- names(selected_columns)
remaining <- setdiff(current_names, "sample_id")
new_order <- c(remaining[1], "sample_id", remaining[-1])
setcolorder(selected_columns, new_order)

selected_columns[,4][is.na(selected_columns[,4]) | selected_columns[,4] == "'--"] <- 0

selected_columns[,3] <- as.numeric(unlist(selected_columns[,3]))

# Convert both columns to numeric safely
selected_columns$demographic.days_to_death <- as.numeric(as.character(unlist(selected_columns$demographic.days_to_death)))
selected_columns$demographic.days_to_birth <- as.numeric(as.character(unlist(selected_columns$demographic.days_to_birth)))

# Check if there are NAs introduced, handle them if needed
selected_columns$demographic.days_to_death[is.na(selected_columns$demographic.days_to_death)] <- 0
selected_columns$demographic.days_to_birth[is.na(selected_columns$demographic.days_to_birth)] <- 0

# Now the arithmetic operation will work
selected_columns$days_difference <- selected_columns$demographic.days_to_death - selected_columns$demographic.days_to_birth

selected_columns <- selected_columns[,-c(3,4)]

selected_columns <- selected_columns[, c(1, 2, 4, 3)]

selected_columns[,4] <- ifelse(selected_columns[,4] == "Alive", 1, 0)

write.table(selected_columns, file="PanCancer_clinical_table_all190708.txt", 
            sep="\t", 
            row.names=FALSE, 
            quote=FALSE)