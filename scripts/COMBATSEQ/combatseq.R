## Install packages as needed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="http://cran.us.r-project.org")
options(repos = BiocManager::repositories())
.bioc_packages <- c("sva", "readr")
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  BiocManager::install(version = "3.18")
  BiocManager::install(.bioc_packages[!.inst], ask = FALSE)
}

library(sva)
library(readr)

# Read environment variables set in the bash script
in_dir <- Sys.getenv("IN_DIR")
out_dir <- Sys.getenv("OUT_DIR")
count_file <- Sys.getenv("COUNT_FILE")



# Read the combined counts file
counts_data <- read_csv(count_file)


# Extract column names excluding 'gene_name'
sample_names <- colnames(counts_data)[-1]  # Assuming the first column is 'gene_name'

# Construct the group variable from the sample names
group <- sample_names

# Extract the biological group information from the sample names
# Assuming the format 'cellType_treatment_batch'
group <- sub("(.*)_[0-9]$", "\\1", sample_names) # This will transform "cDC_134_3" into "cDC_134", "cDC_174_4" into "cDC_174", etc.


# Extract the batch information from the sample names
# Assuming the last character of each sample name is the batch number
batch <- as.numeric(sub(".*_(\\d)$", "\\1", sample_names))

# Check the batch and group variables
# print(head(batch))
# print(head(group))

# Convert counts to matrix if it's not already
counts_matrix <- as.matrix(counts_data[ , -1])  # Assuming the first column is gene identifiers
gene_names <- counts_data$gene_name

# Run ComBat-Seq
combat_adj_data <- ComBat_seq(counts=counts_matrix, batch=batch, group=group)
combat_adj_data_df <- as.data.frame(combat_adj_data)
combat_adj_data_df <- cbind(gene_name = gene_names, combat_adj_data_df)

# Save the adjusted counts
write_csv(combat_adj_data_df, file.path(out_dir, "adjusted_gene_counts.csv"))

print("ComBat-Seq correction completed and output saved.")
