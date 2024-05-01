#!/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/R-4.3.3/bin/Rscript

## Install packages as needed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="http://cran.us.r-project.org")
options(repos = BiocManager::repositories())
.bioc_packages <- c("dplyr", "rtracklayer")
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  BiocManager::install(version = "3.18")
  BiocManager::install(.bioc_packages[!.inst], ask = FALSE)
}

# Load packages
library(dplyr)
library(rtracklayer)

# Define the base directory where data folders are located and the output directory
base_dir <- "/scratch/bioconsult/Ethan_James/BINF_Consulting_2024"
out_dir <- file.path(base_dir, "ongoing/viral_load_quant/PRRSV2")

# Load metadata
metadata <- read.csv(file.path(base_dir, "metadata/PRRSV2_metadata.csv"), header=TRUE)

# Function to calculate TPM
calculate_tpm <- function(count_data, gene_lengths, viral_gene_ids) {
  # Merge counts with gene lengths on Gene IDs
  data_with_lengths <- merge(count_data, gene_lengths, by="Gene")
  
  # Convert counts to RPK
  rpk <- data_with_lengths$Count / (data_with_lengths$Length / 1000)  # Length in kilobases

  # Calculate total RPK to use as denominator in TPM calculation
  total_rpk <- sum(rpk)

  # Calculate TPM for each gene
  tpm <- rpk / total_rpk * 1e6
  tpm_data <- data.frame(Sample=data_with_lengths$Sample, Gene=data_with_lengths$Gene, TPM=tpm)
  filtered_tpm_data <- subset(tpm_data, Gene %in% viral_gene_ids)
  return(filtered_tpm_data)
}

# Function to extract gene lengths from GTF files
get_gene_lengths <- function(gtf_path) {
  gtf <- import(gtf_path)
  genes <- subset(gtf, type == "gene")
  lengths <- width(genes)
  gene_ids <- mcols(genes)$gene_id
  gene_lengths <- data.frame(Gene=gene_ids, Length=lengths)
  row.names(gene_lengths) <- gene_ids
  return(gene_lengths)
}

get_viral_gene_ids <- function(viral_gtf_path) {
    gtf <- import(viral_gtf_path)
    genes <- subset(gtf, type == "gene")
    gene_ids <- mcols(genes)$gene_id
    return(gene_ids)
}

# Initialize a list to store TPM data
all_tpm_data <- list()

filtered_metadata <- metadata[grepl("_Mock_|_134_", metadata$Sample), ]

# Process each sample
for (i in 1:nrow(filtered_metadata)) {
  #strain_name <- ifelse(grepl("134", metadata$Sample[i]), "134", "174")
  strain_name <- ifelse(grepl("Mock", filtered_metadata$Sample[i]), "Mock", "134")
  sample_name <- filtered_metadata$Sample[i]

  combined_gtf_file <- "combined_pig_NC134.gtf"
  viral_gtf_file <- "NC134_annot.gtf"

  gtf_path <- sprintf("%s/raw_data/ref_genomes/PRRSV2/NC134/%s", base_dir, combined_gtf_file)
  viral_gtf_path <- sprintf("%s/raw_data/ref_genomes/PRRSV2/NC134/%s", base_dir, viral_gtf_file)
  
  # Read gene lengths from GTF file
  gene_lengths <- get_gene_lengths(gtf_path)

  # Get viral gene IDs
  viral_gene_ids <- get_viral_gene_ids(viral_gtf_path)
  
  # Construct path to the ReadsPerGene.out.tab file
  file_path <- sprintf("%s/ongoing/mapping/PRRSV2_134_QUANT/PRRSV2_134_QUANT/%s/PRRSV2_134_QUANT_%s_ReadsPerGene.out.tab",
                       base_dir, sample_name, sample_name)
  
  # Read count data
  if (file.exists(file_path)) {
    count_data <- read.delim(file_path, header=FALSE, sep="\t")
    count_data <- count_data[-(1:4), 1:2]
    count_data$Sample <- sample_name 
    # The second column contains the counts
    colnames(count_data) <- c("Gene", "Count", "Sample")
    
    # Calculate TPM
    tpm_data <- calculate_tpm(count_data, gene_lengths, viral_gene_ids)

    print(tpm_data)
    
    # Store in list with metadata
    all_tpm_data[[sample_name]] <- tpm_data
  } else {
    warning(sprintf("File not found for sample %s at path %s", sample_name, file_path))
  }
}

# Combine all TPM data into a single data frame
tpm_results <- do.call(rbind, all_tpm_data)

# Write to file
write.table(tpm_results, file=file.path(out_dir, "viral_load_tpm_results_NC134.tsv"), quote=FALSE, sep="\t", row.names=FALSE)

warnings()

print("Done.")
