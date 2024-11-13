# ==========================================
# Single-Cell RNA-Seq Analysis - LungEndoLPS-Seq
# ==========================================
# Script for DEGs for each cluster 
# LPS vs PBS
# =============== Libraries ================
library(readr)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(Seurat)
library(heatmaply)
library(stringr)
library(plotly)
library(htmltools)
library(data.table)
library(openxlsx)
library(scales)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(RColorBrewer)

# Set working directory for DE analysis
setwd("C:/Users/domin/klaudia_r/DE_analysis_scrnaseq")

# Load integrated Seurat object (assuming this was previously saved)
# intg <- readRDS(file = "path/to/your/file.rds")

# Extract raw counts and metadata from Seurat object to create a SingleCellExperiment object
counts <- intg@assays$RNA@layers$counts 
metadata <- intg@meta.data

# Clean metadata: remove underscores from sample_ID and create combined sample names
metadata$Sample_ID <- gsub("_", "", metadata$Sample_ID)
metadata$combined_name <- paste(metadata$Sample_ID, metadata$sample, sep = "_")

# Create SingleCellExperiment object from Seurat object
sce <- as.SingleCellExperiment(intg)

# Check the assays and counts matrix
dim(counts(sce))
head(counts(sce)[1:6, 1:6])

# Explore cellular metadata
dim(colData(sce))
head(colData(sce))

# Get cluster names and count total clusters
cluster_names <- levels(colData(sce)$ClusterNames)
cat("Total number of clusters:", length(cluster_names), "\n")

# Modify metadata columns for better readability and consistency
colData(sce)$Sample_ID <- gsub("_", "", colData(sce)$Sample_ID)
colData(sce)$combined_name <- paste(colData(sce)$Sample_ID, colData(sce)$sample, sep = "")
sample_names <- unique(colData(sce)$combined_name)
cat("Total number of samples:", length(sample_names), "\n")

# Aggregate counts by cluster and sample group
groups <- colData(sce)[, c("ClusterNames", "combined_name")]
aggr_counts <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum")

# Transpose aggregated matrix for samples as columns and genes as rows
aggr_counts <- t(aggr_counts)

# Inspect structure of column names and split them
tstrsplit(colnames(aggr_counts), "_") %>% str()

# Loop to extract counts for each cluster and store in a list
counts_ls <- list()
for (i in 1:length(cluster_names)) {
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
}

# Check structure of the counts list
str(counts_ls)

# Extract and clean metadata
metadata <- colData(sce) %>% 
  as.data.frame() %>%
  select(sample, Sample_ID, combined_name) %>%
  distinct()

# Ensure metadata row names match the sample names
rownames(metadata) <- metadata$combined_name

# Create a table to count the number of cells per sample and cluster
cell_count_table <- table(colData(sce)$combined_name, colData(sce)$ClusterNames)
cat("Cell counts table (first 6 rows, 7 clusters):\n")
print(cell_count_table[1:6, 1:7])

# Initialize list for metadata specific to each cluster
metadata_ls <- list()

# Loop to match counts with metadata for each cluster
for (i in 1:length(counts_ls)) {
  df <- data.frame(cluster_combined_name = colnames(counts_ls[[i]]))
  df$ClusterNames <- tstrsplit(df$cluster_combined_name, "_")[[1]]
  df$combined_name <- tstrsplit(df$cluster_combined_name, "_")[[2]]
  
  # Get cell count info from global table
  cell_counts <- cell_count_table[, which(colnames(cell_count_table) == unique(df$ClusterNames))]
  cell_counts <- cell_counts[cell_counts > 0] # Exclude zero counts
  
  # Order cell counts and match with sample names
  sample_order <- match(df$combined_name, names(cell_counts))
  df$cell_count <- cell_counts[sample_order]
  
  # Merge with metadata
  df <- plyr::join(df, metadata, by = intersect(names(df), names(metadata)))
  rownames(df) <- df$cluster_combined_name
  
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$ClusterNames)
}

# Save metadata list as RDS file
saveRDS(metadata_ls, file = "metadata_ls.rds")

# Select cell type for differential expression (DE) analysis
selected_cluster <- "LPS induced2"  # Specify cluster of interest
cluster_counts <- counts_ls[[which(names(counts_ls) == selected_cluster)]]
cluster_metadata <- metadata_ls[[which(names(metadata_ls) == selected_cluster)]]

# Set the reference level for DE comparison
cluster_metadata$sample <- factor(cluster_metadata$sample, levels = c("PBS", "LPS"))

# Create DESeq2 object for DE analysis
dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata, design = ~ sample)

# Transform counts for data visualization using rlog
rld <- rlog(dds, blind = TRUE)

# Plot PCA to visualize the sample separation
DESeq2::plotPCA(rld, ntop = 500, intgroup = "sample")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count")

# Compute correlation matrix for rlog-transformed data
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap of correlations
pheatmap(rld_cor, annotation = cluster_metadata[, c("sample"), drop = FALSE])

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Check the results for the DE analysis
resultsNames(dds)

# Generate DE results object
res <- results(dds, name = "sample_LPS_vs_PBS", alpha = 0.05)

# Shrink log2 fold changes using apeglm method for more accurate estimation
res <- lfcShrink(dds, coef = "sample_LPS_vs_PBS")

# Convert results to a tibble for easier handling with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results
head(res_tbl)

# Write DE results to a CSV file
write.csv(res_tbl,
          paste0("Results/", unique(cluster_metadata$ClusterNames), "_", 
                 levels(cluster_metadata$sample)[2], "_vs_", levels(cluster_metadata$sample)[1], "_all_genes.csv"),
          quote = FALSE, row.names = FALSE)

