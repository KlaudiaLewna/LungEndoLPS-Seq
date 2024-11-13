# ==========================================
# Single-Cell RNA-Seq Analysis - LungEndoLPS-Seq
# ==========================================
# Script for GO term enrichment analysis
#DEGS chosen cluster vs gCap with downSample

# =============== Libraries ================
library(dplyr)
library(SeuratObject)
library(Seurat)
library(ggraph)
library(patchwork)
library(grid)
library(gridExtra)
library(DropletUtils)
library(fastmap)
library(SeuratData)
library(SeuratWrappers)
library(BiocManager)
library(multtest)
library(ggpubr)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)
library(caret)
library(clusterProfiler)

# Set options for memory usage
options(future.globals.maxSize = 9912896000)

# Load Seurat object with integrated data
intg <- readRDS("C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\names_res.0.5_6samples_named_rpca.rds")

# Visualize UMAP of the integrated data with cluster labels
DimPlot(intg, reduction = "umap.integrated", label = TRUE)

# Downsampling the data
# Extract cell identities and create a metadata frame
cell_ids <- Cells(intg)
cluster_labels <- Idents(intg)
cell_data <- data.frame(Cell = cell_ids, Cluster = as.factor(cluster_labels))

# Apply downsampling to balance the clusters
downsampled_data <- downSample(x = cell_data, 
                               y = cell_data$Cluster, 
                               list = FALSE, 
                               yname = "Cluster")

# Subset the Seurat object using the downsampled cells
downsampled_cells <- downsampled_data$Cell
intg_subset <- subset(intg, cells = downsampled_cells)

# Differential expression analysis between Cluster 1 and Cluster 2
markers_cluster1_vs_cluster2 <- FindMarkers(intg_subset, ident.1 = "LPS induced1", ident.2 = "gCap", 
                                            min.pct = 0.25, logfc.threshold = 0.25)

# View top results from the differential expression analysis
head(markers_cluster1_vs_cluster2)

# Filter significant genes (adjusted p-value < 0.05)
significant_genes <- markers_cluster1_vs_cluster2 %>%
  filter(p_val_adj < 0.05)

# Extract the gene names for upregulated and downregulated genes
upregulated_genes <- rownames(significant_genes[significant_genes$avg_log2FC > 0, ])
downregulated_genes <- rownames(significant_genes[significant_genes$avg_log2FC < 0, ])

# Perform GO enrichment analysis for Biological Processes (BP) on upregulated genes
go_enrichment_up <- enrichGO(gene = upregulated_genes, 
                             OrgDb = org.Mm.eg.db,  # Change to org.Hs.eg.db for human data
                             keyType = "SYMBOL", 
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05)

# View the GO enrichment results
head(go_enrichment_up)

# Generate a bar plot for the GO terms
barplot(go_enrichment_up, showCategory = 20)

# Perform GO enrichment analysis again with clusterProfiler
ego <- enrichGO(gene = upregulated_genes,
                OrgDb = org.Mm.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

# Filter for transcription-related GO terms
transcription_terms <- ego@result[grep("transcription|RNA|DNA binding", ego@result$Description, ignore.case = TRUE),]

# View or save transcription-related GO terms
print(transcription_terms)

# Create an enrichment map (network of GO terms)
# emapplot(go_enrichment_up)

# Heatmap of Differential Expression
# Select top 50 differentially expressed genes based on adjusted p-value
top50_genes <- markers_cluster1_vs_cluster2 %>%
  arrange(p_val_adj) %>%
  head(50) %>%
  rownames()

# Subset Seurat object to include only Cluster 1 and Cluster 2
intg_clusters <- subset(intg_subset, idents = c("LPS induced1", "gCap"))

# Calculate the average expression for each gene in the selected clusters
avg_expression <- AverageExpression(intg_clusters, features = top50_genes)$SCT

# Log-normalize the average expression values
log_normalized_expression <- log1p(avg_expression)

# Plot a heatmap of log-normalized values
pheatmap(as.matrix(log_normalized_expression), 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Log-Normalized Average Expression of Top Differentially Expressed Genes Between Cluster 1 and Cluster 2")
