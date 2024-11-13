# ==========================================
# Single-Cell RNA-Seq Analysis - LungEndoLPS-Seq
# ==========================================
# Script for doublet identification,.
# =============== Libraries ================
# Required Libraries
library(dplyr)
library(sp)
library(SeuratObject)
library(Seurat)
library(sctransform)
library(clustree)
library(DoubletFinder)
library(stringr)
library(ggraph)
library(ggplot2)
library(patchwork)
library(SeuratWrappers)

# Set maximum memory size for future operations
options(future.globals.maxSize = 891289600099)

# Load Seurat Object
LPS <- readRDS(file = "C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2020_024\\LPS_2020_024_Redy")

# Steps for each sample were the same
# ----------------------------------------------
# Step 1: DoubletFinder pK Identification
# ----------------------------------------------
# Perform pK sweep to identify optimal pK value for doublet detection
sweep.res.list_LPS <- paramSweep(LPS, PCs = 1:30, sct = TRUE)
sweep.stats_LPS <- summarizeSweep(sweep.res.list_LPS, GT = FALSE)
bcmvn_LPS <- find.pK(sweep.stats_LPS)

# Plot pK vs BCmetric to select optimal pK
ggplot(bcmvn_LPS, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

# Select the pK with the highest BCmetric value
pK <- bcmvn_LPS %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) %>%
  as.numeric()

# ----------------------------------------------
# Step 2: Check Sample_ID and Cell Counts
# ----------------------------------------------
# Check if "Sample_ID" exists and count cells per sample
if ("Sample_ID" %in% colnames(LPS@meta.data)) {
  cell_counts <- LPS@meta.data %>%
    group_by(Sample_ID) %>%
    summarize(Cell_Count = n())
  print(cell_counts)
} else {
  print("Sample_ID column not found in metadata.")
}

# ----------------------------------------------
# Step 3: Homotypic Doublet Proportion Estimate
# ----------------------------------------------
if ("seurat_clusters" %in% colnames(LPS@meta.data)) {
  homotypic.prop <- modelHomotypic(LPS$seurat_clusters)
  nExp_poi <- round(0.056 * nrow(LPS@meta.data))  # Adjust according to your doublet rate
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder with adjusted parameters
  LPS <- doubletFinder(LPS, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
  
  # Plot UMAP with DoubletFinder classifications
  df_col_name <- paste0("DF.classifications_0.25_", pK, "_", nExp_poi.adj)
  if (df_col_name %in% colnames(LPS@meta.data)) {
    DimPlot(LPS, reduction = 'umap', group.by = df_col_name)
  } else {
    print(paste("Column", df_col_name, "not found in metadata."))
  }
} else {
  print("seurat_clusters column not found in metadata.")
}

# ----------------------------------------------
# Step 4: Count Singlets and Doublets
# ----------------------------------------------
# Check classification for singlets and doublets
table(LPS@meta.data$DF.classifications_0.25_0.11_286)

# Count doublets per cluster
doublet_per_cluster <- LPS@meta.data %>%
  group_by(SCT_snn_res.0.5, DF.classifications_0.25_0.15_744) %>%
  summarize(doublet_count = n()) %>%
  filter(DF.classifications_0.25_0.15_744 == "Doublet")
print(doublet_per_cluster)

# ----------------------------------------------
# Step 5: Filter Out Doublets
# ----------------------------------------------
# Function to filter out doublets
filter_doublets <- function(seurat_obj, df_col_name) {
  if (!df_col_name %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", df_col_name, "not found in metadata."))
  }
  singlets <- seurat_obj[, seurat_obj@meta.data[[df_col_name]] == "Singlet"]
  return(singlets)
}

# Specify the classification column and filter doublets
df_col_name <- "DF.classifications_0.25_0.11_286"
LPS_singlets <- filter_doublets(LPS, df_col_name)

# Print number of cells before and after filtering
num_cells_before <- ncol(LPS)
num_cells_after <- ncol(LPS_singlets)
cat("Number of cells before filtering:", num_cells_before, "\n")
cat("Number of cells after filtering:", num_cells_after, "\n")

# Save full dataset with doublet classifications
saveRDS(LPS, file = "C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2020_024\\LPS_20_024_doublet_finder.rds")

# Save filtered data
saveRDS(LPS_singlets, file = "C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2020_024\\LPS_20_024_removed_doublets.rds")


# ----------------------------------------------
# Step 6: Filter Clusters Based on Doublet Proportions
# ----------------------------------------------
# Read filtered Seurat object and visualize UMAP
LPS <- readRDS("C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2024_024\\LPS_24_024_doublet_finder.rds")
LPS_singlets <- readRDS("C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2024_024\\LPS_24_024_removed_doublets.rds")
DimPlot(LPS_singlets, reduction = 'umap', group.by = "SCT_snn_res.0.5")

# Remove specific clusters (e.g., cluster 12)
filtered_LPS <- subset(LPS_singlets, subset = SCT_snn_res.0.5 == "8", invert = TRUE)

# Check number of cells per cluster after filtering
cell_counts_per_cluster <- table(Idents(filtered_LPS))
print(cell_counts_per_cluster)

# Save filtered Seurat object
saveRDS(filtered_LPS, file = "C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2024_026\\LPS_24_026_removed_doublets_badclust.rds")

# ----------------------------------------------
# Step 7: Marker Gene Analysis
# ----------------------------------------------
# Define and plot top markers for different clusters and cell types
top_markers <- c()
LPS.markers <- FindAllMarkers(LPS, only.pos = TRUE)
for (i in 0:(length(unique(Idents(LPS))) - 1)) {
  top_markers <- c(top_markers, head(subset(LPS.markers, LPS.markers$cluster == i)$gene, 5))
}
top_markers <- top_markers[-c(31, 83)]  # Remove unwanted markers
DotPlot(LPS, features = top_markers) + RotatedAxis() + labs(title = "Top markers for clusters")

# Define marker lists for specific cell types
EC_immune_Mlana_markers <- list("EC" = c("Pecam1", "Cdh5", "Tie1", "Tie2", "Kdr"),
                                "Immune" = c("Ncam1", "Ptprc", "Cd68", "Fcgr3", "Csf3r", "Cd8a", "Cd4", "Cd19"))
Niina_markers <- list("Artery" = c("Adgrg6", "Mecom", "Bmx", "Thsd7a", "Atp13a3", "Ltbp4", "Plac8", "Mgst1"),
                      "gCap" = c("Kit", "Sema3g", "Lpl", "Hey1", "Hmcn1", "Glp1r", "Cadm1", "Epb41"))

# Plot specific marker gene expression
DotPlot(LPS, features = Niina_markers) + RotatedAxis() + labs(title = "8 markers")
DotPlot(LPS, features = EC_immune_Mlana_markers) + RotatedAxis() + labs(title = "Markers for EC and immune cells")

# Save analysis results
saveRDS(LPS, file = "LPS_final_analysis.rds")

