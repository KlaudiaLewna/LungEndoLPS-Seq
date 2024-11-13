# ==========================================
# Single-Cell RNA-Seq Analysis - LungEndoLPS-Seq
# ==========================================
# Script for integration of all samples
# elimination of batch effect
# =============== Libraries ================

# Load necessary libraries
library(dplyr)
library(Seurat)
library(sctransform)
library(SeuratWrappers)
library(clustree)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)
library(SoupX)
library(DropletUtils)

# Set options
options(future.globals.maxSize = 89128960463630)

# Function for calculating median SDs
get_median_sds_as_v <- function(x, two_sd = TRUE) {
  # Function to compute median SDs
  # Add custom logic here if needed
}

# Loading data
LPS_20_024 <- readRDS(file = "C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2020_024\\LPS_20_024_removed_doublets.rds")
LPS_20_025 <- readRDS(file = "C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2020_025\\LPS_20_025_removed_doublets.rds")
LPS_24_024 <- readRDS(file = "C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2024_024\\LPS_24_024_removed_doublets_badclust.rds")
LPS_24_025 <- readRDS(file = "C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2024_025\\LPS_24_025_removed_doublets_badclust.rds")
LPS_24_026 <- readRDS(file = "C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2024_026\\LPS_24_026_removed_doublets_badclust.rds")
LPS_24_027 <- readRDS(file = "C:\\Users\\domin\\klaudia_r\\Data_BetterQuality\\2024_027\\LPS_24_027_removed_doublets_badclust.rds")


# Merge data
LPS <- merge(x = LPS_20_024, y = list(LPS_20_025, LPS_24_024, LPS_24_025, LPS_24_026, LPS_24_027))

# Save merged data for future use
saveRDS(LPS, "C:\\Users\\domin\\klaudia_r\\Cancer_controls\\merged_6controls.rds")
LPS <- readRDS("C:\\Users\\domin\\klaudia_r\\Cancer_controls\\merged_6controls.rds")

# Add 'sample' column
LPS$sample <- ifelse(grepl("LPS", rownames(LPS)), "LPS", "PBS")

# Perform SCTransform
LPS <- SCTransform(LPS, vst.flavor = "v2")

# PCA analysis
LPS <- RunPCA(LPS)

# Identify cell clusters (no integration)
obj <- LPS
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, graph.name = "unintegrated")
obj <- FindClusters(obj, resolution = 0.5, algorithm = 1, graph.name = "unintegrated")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated")
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("Sample_ID", "unintegrated_clusters"))

# Data integration using RPCA
obj2 <- IntegrateLayers(
  object = LPS,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  normalization.method = "SCT",
  verbose = TRUE
)

# Find clusters after integration
obj2 <- FindNeighbors(obj2, reduction = "integrated.rpca", dims = 1:30)
for (res in c(0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4)) {
  obj2 <- FindClusters(obj2, resolution = res, algorithm = 1, cluster.name = paste("rpca_clusters_res.", as.character(res), sep = ""), verbose = FALSE)
}

# Plot clustering tree
clustree(obj2@meta.data, prefix = "rpca_clusters_res.")

# Run UMAP for visualization
obj2 <- RunUMAP(obj2, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.integrated")
DimPlot(obj2, reduction = "umap.integrated", group.by = c("Sample_ID", "rpca_clusters_res.0.5"), label = FALSE)

# Save integrated object with RPCA
obj2$rpca_clusters <- obj2$rpca_clusters_res.0.5
Idents(obj2) <- obj2$rpca_clusters_res.0.5
saveRDS(obj2, "C:\\Users\\domin\\klaudia_r\\Cancer_controls\\integrated_6ctrl_samples_rpca_better.rds")

# Harmony integration
obj4 <- IntegrateLayers(
  object = LPS,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony",
  normalization.method = "SCT",
  verbose = FALSE
)

# Find clusters after Harmony integration
obj4 <- FindNeighbors(obj4, reduction = "integrated.harmony", dims = 1:30)
for (res in c(0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4)) {
  obj4 <- FindClusters(obj4, resolution = res, algorithm = 1, cluster.name = paste("harmony_clusters_res.", as.character(res), sep = ""), verbose = FALSE)
}

# Plot clustering tree for Harmony
clustree(obj4@meta.data, prefix = "harmony_clusters_res.")

# Run UMAP for Harmony integration
obj4 <- RunUMAP(obj4, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.integrated")
DimPlot(obj4, reduction = "umap.integrated", group.by = c("Sample_ID", "harmony_clusters_res.0.5"))

# Save integrated object with Harmony
obj4$harmony_clusters <- obj4$harmony_clusters_res.0.5
Idents(obj4) <- obj4$harmony_clusters_res.0.5
saveRDS(obj4, "C:\\Users\\domin\\klaudia_r\\Cancer_controls\\integrated_6ctrl_harmony.rds")
