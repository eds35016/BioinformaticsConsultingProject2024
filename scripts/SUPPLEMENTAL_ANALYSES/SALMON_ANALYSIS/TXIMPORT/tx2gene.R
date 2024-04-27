# Set CRAN repository
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Check and install BiocManager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("tximport")


# Check and install readr if necessary
if (!requireNamespace("readr", quietly = TRUE))
    install.packages("readr")


library(tximport)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
quant_file <- args[1]
output_dir <- args[2]
mapping_file <- args[3]

# Read the transcripts-to-gene mapping file
tx2gene <- read.delim(mapping_file, header=FALSE, col.names=c("transcript", "gene"))


# Import transcript-level counts and aggregate to gene-level
txi <- tximport(quant_file, type = "salmon", tx2gene = tx2gene)

# Access the gene-level counts
gene_counts <- txi$counts


# Write the gene-level count matrix to a csv
write.csv(as.data.frame(gene_counts), file=file.path(output_dir, "gene_counts_matrix.csv"))


