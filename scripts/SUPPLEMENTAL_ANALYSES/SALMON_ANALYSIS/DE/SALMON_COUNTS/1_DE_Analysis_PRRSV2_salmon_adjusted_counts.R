#!/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/R-4.3.3/bin/Rscript
## Install packages as needed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="http://cran.us.r-project.org")
options(repos = BiocManager::repositories())
.bioc_packages <- c("edgeR", "limma", "ComplexHeatmap", "circlize", "svglite", "dplyr", "stringr", "ggplot2")
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  BiocManager::install(version = "3.18")
  BiocManager::install(.bioc_packages[!.inst], ask = FALSE)
}

library(edgeR)
library(ComplexHeatmap)
library(ggplot2)
library(limma)
library(circlize)
library(svglite)
library(dplyr)
library(stringr)

metadata_dir <- "/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/metadata"
count_dir <- "/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/PRRSV2/PRRSV2"
out_dir <- "/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/DE/salmon_counts/PRRSV2_adjusted_postQC"

metadata <- read.csv(file.path(metadata_dir, "PRRSV2_metadata.csv"))
metadata <- metadata[!metadata$Sample %in% c("cDC_134_3", "moDC_Mock_7"), ] # Eliminate QC removed samples from metadata
count_matrix <- read.csv(file.path(count_dir, "adjusted_gene_counts_salmon_postQC.csv"), row.names = 1)
metadata$Animal <- as.numeric(gsub("\\D", "", metadata$Sample))
metadata$CellType <- str_extract(metadata$Sample, "^[^_]*")
metadata$Strain <- str_extract(metadata$Sample, "(?<=_)[^_]*(?=_)")

de_results <- list()

for (CellType in unique(metadata$CellType)) {

  # Filter metadata_series to only include the current cell type
  metadata_CT <- metadata[metadata$CellType == CellType, ]
  print(metadata_CT)

  # Filter count matrix to only include the current cell type
  count_matrix_CT <- count_matrix[, colnames(count_matrix) %in% metadata_CT$Sample]

  # Keep column ordering the same as metadata
  count_matrix_CT_reordered <- count_matrix_CT[, metadata_CT$Sample]
  count_matrix_CT <- count_matrix_CT_reordered

  metadata_CT$Treatment <- factor(metadata_CT$Treatment)
  metadata_CT$Treatment <- relevel(metadata_CT$Treatment, ref = paste0(CellType, "_Mock"))

  design_matrix <- model.matrix(~0+Treatment, data=metadata_CT)
  colnames(design_matrix) <- levels(metadata_CT$Treatment)
  colnames(design_matrix) <- make.names(colnames(design_matrix))

  print(design_matrix)

  keep <- filterByExpr(y = count_matrix_CT, design = design_matrix)
  count_matrix_CT <- count_matrix_CT[keep, ]
  y <- DGEList(counts = count_matrix_CT)
  y <- calcNormFactors(y)

  ###########################################
  # Limma Trend Approach

  #logCPM <- cpm(y, log = TRUE, prior.count = 3)

  #fit <- lmFit(logCPM, design_matrix)
  #fit <- eBayes(fit, trend=TRUE)

  ###########################################

  # Limma Voom Approach

  v <- voom(y, design_matrix, plot=TRUE)

  fit <- lmFit(v, design_matrix)
  fit <- eBayes(fit)
  topTable(fit, coef=ncol(design_matrix))

  ###########################################

  # Initialize an empty list to store named contrast expressions
  namedContrastsList <- list()

  # Loop through treatments to create and name each contrast
  for(Treatment in metadata_CT$Treatment) {
    strain <- unique(metadata_CT$Strain[metadata_CT$Treatment == Treatment])
    if(strain != "Mock") {
      # Define the contrast expression
      contrastExpr <- paste0(Treatment, " - ", CellType, "_Mock")
      contrastName <- make.names(paste0(Treatment, ".vs.Mock"), unique = TRUE)
      
      # Assign the contrast expression to the named list
      namedContrastsList[[contrastName]] <- contrastExpr
    }
  }

  # Add magnitude comparison
  contrastExpr <- paste0("(", namedContrastsList[2], ")-(", namedContrastsList[1], ")")
  contrastName <- make.names(paste0(CellType, ".134vMock.vs.174vMock"), unique = TRUE)
  namedContrastsList[[contrastName]] <- contrastExpr
  print(paste0(CellType, "- ", contrastName, " ", contrastExpr))

  # Ensure 'levels' argument is correctly passed
  namedContrastsList[["levels"]] = design_matrix

  # Use 'do.call' to dynamically call 'makeContrasts' with the list of named contrasts
  contrasts <- do.call(makeContrasts, namedContrastsList)

  fit2 <- contrasts.fit(fit, contrasts)

  # Trend Approach
  #fit2 <- eBayes(fit2, trend=TRUE)

  # Voom Approach
  fit2 <- eBayes(fit2)

  print(head(fit2))

  de.m <- NULL;
  lfc.ct <- log2(2.5);
  p.ct <- 0.05;
  for (i in 1:ncol(contrasts)) {
    coef.name <- colnames(contrasts)[i];
    t <- topTable(fit2, coef=coef.name, number=nrow(fit2), sort="none");
    t <- t[,c("logFC","AveExpr", "P.Value", "adj.P.Val")];
    flag <- abs(t[,"logFC"]) > lfc.ct & t[,"P.Value"] < p.ct; # raw p-values
    t <- data.frame(t, flag, stringsAsFactors=F);

    colnames(t) <- paste(coef.name, colnames(t), sep="."); # label the columns
    if (is.null(de.m)) { de.m <- t;
    } else {
      if (!identical(rownames(de.m), rownames(t))) stop("error");
      de.m <- cbind(de.m, t);
    }
  }

  # MA Plots
  # MA Plot generation function
  plotMA_contrast <- function(fit_obj, coef_index, lfc_threshold = log2(2.5), pval_threshold = 0.05) {
    res <- topTable(fit_obj, coef=coef_index, number=Inf, sort.by="none")
    
    p <- ggplot(res, aes(x=AveExpr, y=logFC)) + 
      geom_point(aes(color=adj.P.Val < pval_threshold & abs(logFC) > lfc_threshold), alpha=0.4) +
      geom_smooth(method="loess", se=FALSE, color="black", linetype="dotted") + 
      scale_color_manual(values = c("grey", "red")) +
      geom_vline(xintercept=0, linetype="dashed", color="blue") + 
      geom_hline(yintercept=c(-lfc_threshold, lfc_threshold), linetype="dashed", color="red") + 
      theme_bw() +
      labs(title=paste("MA Plot:", coef_index),
          x="Average Expression", y="Log Fold Change")

    return(p)
  }

  # Generate and save MA plots for each contrast
  contrast_names <- colnames(contrasts)
  for (contrast_name in contrast_names) {
    svg(filename=file.path(out_dir, paste0("MA_plot_", contrast_name, ".svg")), width=15, height=5)
    print(plotMA_contrast(fit2, contrast_name, lfc.ct, p.ct))
    dev.off()
  }

  # Identify which columns end in '.flag'
  flag_cols <- grep("\\.flag$", colnames(de.m), value = TRUE)

  # Save differential expression matrix de.m to tsv file
  write.table(de.m, file = file.path(out_dir, paste0(CellType, "_de.m.tsv")), sep = "\t")

  # Subset the de.m data frame to only include significant genes
  sig_genes <- de.m[apply(de.m[, flag_cols], 1, any), ]

  print("Total number of significantly differentially expressed genes:")
  print(nrow(sig_genes))

  # Check if there are any significant genes based on logFC and adjusted P-Value
  if(nrow(sig_genes) == 0) {
    warning("No significant genes found with adjusted p-value < 0.05 and abs(logFC) > log2(2.5) in any of the contrasts")

  }

  # Trend Approach
  #sig_counts <- logCPM[rownames(logCPM) %in% rownames(sig_genes), ]

  # Voom Approach
  sig_counts_full <- v[rownames(v) %in% rownames(sig_genes), ]
  sig_counts <- sig_counts_full$E

  write.table(sig_counts, file = file.path(out_dir, paste0(CellType, "_sig_counts.tsv")), sep = "\t")

  ######################################################

  all_counts <- v$E

  # Create a color mapping
  min_val <- -4
  max_val <- 4
  color_mapping <- colorRamp2(c(min_val, 0, max_val), c("blue", "white", "red"))

  #######################################

  # Extract logFC columns for significant genes
  logFC_matrix <- as.matrix(de.m[rownames(de.m) %in% rownames(sig_genes), grep("logFC", colnames(de.m))])

  # Create a factor that will determine column splits
  column_split_factor <- factor(colnames(contrasts), levels = unique(colnames(contrasts)))

  # Treatment heatmap
  svglite(file.path(out_dir, paste0(CellType, "_Treatment_Heatmap.svg")), width=8, height=14)
  p <- Heatmap(
      logFC_matrix, 
      name = "logFC", 
      show_row_names = FALSE, 
      show_column_names = TRUE,
      cluster_rows = TRUE, 
      cluster_columns = FALSE,
      col=color_mapping,
      column_split = column_split_factor,
      column_title = NULL
  )

  print(p)
  dev.off()

  print(metadata_CT$Treatment)

  # MDS Plot
  pdf(file.path(out_dir, paste0(CellType, "_MDS_Plot.pdf")), height=10, width=10)
  plotMDS(y, col=as.integer(metadata_CT$Treatment))
  legend("topright", legend=levels(metadata_CT$Treatment), col=1:length(levels(metadata_CT$Treatment)), pch=20)
  dev.off()
}