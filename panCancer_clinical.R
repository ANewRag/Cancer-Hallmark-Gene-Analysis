library(data.table)

# Read in the sample and clinical sheets
sample_sheet <- read.delim("gdc_sample_sheet.2025-04-18 (3).tsv", 
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
clinical_sheet <- fread("clinical.tsv")

# Remove duplicated rows based on the second column
removed_dup <- clinical_sheet[!duplicated(clinical_sheet[[2]]), ]

# Subset the clinical sheet to keep the desired columns:
#   Column 10: original case id (for matching)
#   Column 11: additional column to keep
#   Column 17: demographic.days_to_death
#   Column 18: demographic.days_to_birth
#   Column 22: additional column to keep
#   Column 27: additional column to keep
#   Column 29: demographic.vital_status (vital status)
selected_columns <- removed_dup[, .SD, .SDcols = c(10, 11, 17, 18, 22, 27, 29)]

setDT(sample_sheet)

# Add the sample_id column by matching the original case id to sample_sheet
selected_columns[, sample_id := sample_sheet[match(get(names(selected_columns)[1]), Case.ID), Sample.ID]]

# Remove the original case id column (first column) as it is no longer needed
orig_col <- names(selected_columns)[1]
selected_columns[, (orig_col) := NULL]

# Reorder columns so that sample_id comes first
setcolorder(selected_columns, c("sample_id", setdiff(names(selected_columns), "sample_id")))

# Replace NA or "'--" in the demographic.vital_status column with 0
selected_columns[is.na(demographic.vital_status) | demographic.vital_status == "'--", 
                 demographic.vital_status := 0]

# Convert the demographic.days_to_death and demographic.days_to_birth columns to numeric
selected_columns[, demographic.days_to_death := as.numeric(as.character(demographic.days_to_death))]
selected_columns[, demographic.days_to_birth := as.numeric(as.character(demographic.days_to_birth))]
selected_columns[is.na(demographic.days_to_death), demographic.days_to_death := 0]
selected_columns[is.na(demographic.days_to_birth), demographic.days_to_birth := 0]

# Create the days_difference column (days_to_death minus days_to_birth)
selected_columns[, days_difference := demographic.days_to_death - demographic.days_to_birth]

# Remove the original numeric demographic columns as they are no longer needed
selected_columns[, c("demographic.days_to_death", "demographic.days_to_birth") := NULL]

# Rearrange columns so that the final order is: sample_id, days_difference, demographic.vital_status,
# and then any remaining columns
setcolorder(selected_columns, c("sample_id", "days_difference", "demographic.vital_status",
                                setdiff(names(selected_columns), c("sample_id", "days_difference", "demographic.vital_status"))))

# Convert the demographic.vital_status column to binary: if the value is "Alive", set it to 1, otherwise 0
selected_columns[, demographic.vital_status := ifelse(demographic.vital_status == "Alive", 1, 0)]

# Rename the gender and race columns if they exist
if ("Gender(0=female;1=male)" %in% names(selected_columns)) {
  setnames(selected_columns, "Gender(0=female;1=male)", "demographic.gender")
}
if ("Race(1=white;2=asian;3=black/african american)" %in% names(selected_columns)) {
  setnames(selected_columns, "Race(1=white;2=asian;3=black/african american)", "demographic.race")
}

# Convert the demographic.gender and demographic.race columns to their numerical representations.
# For gender: "female" becomes 0, "male" becomes 1.
if ("demographic.gender" %in% names(selected_columns)) {
  selected_columns[, demographic.gender := ifelse(demographic.gender == "female", 0,
                                                  ifelse(demographic.gender == "male", 1, NA))]
}

# For race: "white" becomes 1, "asian" becomes 2, and "black/african american" becomes 3.
if ("demographic.race" %in% names(selected_columns)) {
  selected_columns[, demographic.race := ifelse(demographic.race == "white", 1,
                                                ifelse(demographic.race == "asian", 2,
                                                       ifelse(demographic.race == "black/african american", 3, NA)))]
}

# Write the final table to a file
write.table(selected_columns, file = "PanCancer_clinical_table_all190708.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
