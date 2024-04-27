#!/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/R-4.3.3/bin/Rscript
# Set CRAN repository
options(repos = c(CRAN = "https://cloud.r-project.org/"))


# Load required libraries or install them if necessary
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(readr)
library(tibble)
library(dplyr)

# Set the working directory to where your subdirectories are located
# Replace this with your actual directory path
setwd("/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/PRRSV2/PRRSV2")

# List all subdirectories containing the count files
subdirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)[-1]  # Excludes the current directory '.'

# Initialize an empty list to store data from each file
data_list <- list()

# Initialize an empty dataframe for the combind data
combined_data <- tibble()

for (dir in subdirs) {
  file_path <- file.path(dir, "gene_counts_matrix.csv")
  # Extract sample name from the directory name
  sample_name <- basename(dir)
  
  # Check if the file exists before trying to read it - ALSO remove two problematic samples
  if (file.exists(file_path) && sample_name != "cDC_134_3" && sample_name != "moDC_Mock_7") {
    # Read the counts file, handle the first column as row names
    counts_data <- read_csv(file_path, col_names = c("gene_name", "count"), skip=1,show_col_types = FALSE)
    counts_data
    
    # Rename the count column to match the sample name
    counts_data <- counts_data %>%
      rename(!!sample_name := count)
    counts_data
    # If combined_data is empty, initialize it with the first file's data
    if (nrow(combined_data) == 0) {
      combined_data <- counts_data
    } else {
      # Otherwise, merge this file's data with the combined data by gene name
      combined_data <- combined_data %>%
        full_join(counts_data, by = "gene_name")
    }
  }
}

# Write the final combined data frame to a CSV file
write_csv(combined_data, "combined_gene_counts.csv")

# Print a message to indicate completion
print("Combined gene counts data frame has been saved as 'combined_gene_counts.csv'.")

