library(data.table)

# Read in the sample and clinical sheets
sample_sheet <- read.delim("gdc_sample_sheet.2025-04-06.tsv", 
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
clinical_sheet <- fread("clinical.tsv")

# Remove duplicated rows based on the second column
removed_dup <- clinical_sheet[!duplicated(clinical_sheet[[2]]), ]

# Subset the clinical sheet to keep the desired columns:
#   Column 10: (original case id, to be used only for matching)
#   Column 17: demographic.days_to_death
#   Column 18: demographic.days_to_birth
#   Column 29: status (for example, "Alive")
selected_columns <- removed_dup[, .SD, .SDcols = c(10, 17, 18, 29)]

setDT(sample_sheet)

# Add the sample_id column by matching the original case id to sample_sheet
selected_columns[, sample_id := sample_sheet[match(get(names(selected_columns)[1]), Case.ID), Sample.ID]]

# Remove the original case id column (from clinical sheet) as it is no longer needed
orig_col <- names(selected_columns)[1]
selected_columns[, (orig_col) := NULL]

# Reorder columns so that sample_id comes first;
# the remaining columns come in their original order
setcolorder(selected_columns, c("sample_id", setdiff(names(selected_columns), "sample_id")))

# At this point, assume the remaining columns are:
#   Column 2: demographic.days_to_death
#   Column 3: demographic.days_to_birth
#   Column 4: the status column (e.g. vital status, possibly holding "Alive")
# For the status column, replace NA or the string "'--" with 0
selected_columns[, 4][is.na(selected_columns[, 4]) | selected_columns[, 4] == "'--"] <- 0

# Convert the demographic columns to numeric
selected_columns$demographic.days_to_death <- as.numeric(as.character(selected_columns$demographic.days_to_death))
selected_columns$demographic.days_to_birth <- as.numeric(as.character(selected_columns$demographic.days_to_birth))
selected_columns$demographic.days_to_death[is.na(selected_columns$demographic.days_to_death)] <- 0
selected_columns$demographic.days_to_birth[is.na(selected_columns$demographic.days_to_birth)] <- 0

# Create the days_difference column (days_to_death minus days_to_birth)
selected_columns[, days_difference := demographic.days_to_death - demographic.days_to_birth]

# Remove the original numeric columns that are no longer needed
selected_columns[, c("demographic.days_to_death", "demographic.days_to_birth") := NULL]

# Now, the remaining columns should be: sample_id, status (from clinical_sheet originally column 29),
# and days_difference (just computed).
# Rearrange columns so that the final order is: sample_id, days_difference, status.
setcolorder(selected_columns, c("sample_id", "days_difference", setdiff(names(selected_columns), c("sample_id", "days_difference"))))

# Convert the status column to binary: if the value is "Alive", set it to 1, otherwise 0.
# (Here we assume the status column’s name is the only one not "sample_id" or "days_difference".)
status_col <- setdiff(names(selected_columns), c("sample_id", "days_difference"))
selected_columns[, (status_col) := ifelse(get(status_col) == "Alive", 1, 0)]

# Write the final table to a file
write.table(selected_columns, file = "PanCancer_clinical_table_all190708.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)