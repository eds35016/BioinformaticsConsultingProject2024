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
setwd("/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/PRRSV2/PRRSV2")

# List all subdirectories containing the count files
subdirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)[-1]  # Excludes the current directory '.'

# Initialize an empty list to store data from each file
data_list <- list()

# Initialize an empty dataframe for the combind data
combined_data <- tibble()

for (dir in subdirs) {
  # Construct the file path with the new naming pattern
  sample_name <- basename(dir)
  file_name <- paste0("PRRSV2_", sample_name, "_ReadsPerGene.out.tab")
  file_path <- file.path(dir, file_name)

  # Check if the file exists before trying to read it - ALSO remove two problematic samples
  if (file.exists(file_path) && sample_name != "cDC_134_3" && sample_name != "moDC_Mock_7") {
    # Read the file, skipping the first four summary lines and only the first two columns
    counts_data <- read_tsv(file_path, col_names = c("gene_name", "count"), skip = 4, col_types = cols_only(gene_name = col_character(), count = col_double()), show_col_types = FALSE)
    # Rename the count column to match the sample name
    counts_data <- counts_data %>%
      rename(!!sample_name := count)
    
    # Merge with the combined data
    if (nrow(combined_data) == 0) {
      combined_data <- counts_data
    } else {
      combined_data <- combined_data %>%
        full_join(counts_data, by = "gene_name")
    }
  }
}

# Write the final combined data frame to a CSV file
write_csv(combined_data, "combined_gene_counts.csv")

print("Combined gene counts data frame has been saved as 'combined_gene_counts.csv'.")
