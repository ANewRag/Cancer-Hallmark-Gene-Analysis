# Load required library
library(fs)

flatten_directories <- function(source_folder, target_folder) {
  # Check if the target folder exists, if not, create it
  if (!dir_exists(target_folder)) {
    dir_create(target_folder)
  }
  
  # List all files in source_folder and subfolders
  files <- dir_ls(source_folder, recurse = TRUE, type = "file")
  
  # Loop through each file
  for (file in files) {
    # Get the file name (basename of the file path)
    file_name <- basename(file)
    
    # Define the target path
    target_path <- file.path(target_folder, file_name)
    
    # Handle file name collisions (if the file already exists in the target folder)
    if (file_exists(target_path)) {
      # Extract the name and extension
      file_ext <- tools::file_ext(file_name)
      base_name <- tools::file_path_sans_ext(file_name)
      counter <- 1
      
      # Generate a new unique name if the file already exists
      while (file_exists(target_path)) {
        new_file_name <- paste0(base_name, "_", counter, ".", file_ext)
        target_path <- file.path(target_folder, new_file_name)
        counter <- counter + 1
      }
    }
    
    # Move the file to the target folder
    file_copy(file, target_path)
    message("Moved: ", file, " -> ", target_path)
  }
}

# Example Usage
source_folder <- "/Users/junninghu/Downloads/gdc_download_20250328_001817.786714_2"  # Replace with the path to the source folder containing subfolders
target_folder <- "/Users/junninghu/Downloads/gdc_download_20250328_001817.786714_2"  # Replace with the path to the folder where you want all the files

flatten_directories(source_folder, target_folder)
