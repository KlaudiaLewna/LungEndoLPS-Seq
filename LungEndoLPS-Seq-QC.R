# ==========================================
# Single-Cell RNA-Seq Analysis - LungEndoLPS-Seq
# ==========================================
# Script for loading, processing, filtering, and visualizing scRNA-seq data.
# This code prepares the data for downstream analyses, including quality control,
# filtering, and subsetting based on treatment.

# =============== Libraries ================
# Load required libraries
library(SeuratObject)
library(Seurat)
library(sctransform)
library(clustree)
library(DoubletFinder)
library(ggraph)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)
library(DropletUtils)
library(fastmap)
library(SeuratWrappers)

# Set global options
options(future.globals.maxSize = 891289600)

# Function to calculate threshold values based on median and standard deviation
get_median_sds_as_v <- function(value, two_sd) {
  sds <- c(median(value) - 5 * sd(value),
           median(value) - 4 * sd(value),
           median(value) - 3 * sd(value),
           median(value),
           median(value) + 3 * sd(value),
           median(value) + 4 * sd(value),
           median(value) + 5 * sd(value))
  
  if (two_sd) {
    sds <- c(sds[1:3], median(value) - 2 * sd(value), sds[4], median(value) + 2 * sd(value), sds[5:7])
  }
  return(sds)
}

# Function to display Seurat cluster information in a table
tabulate_seurat_clusters_info <- function(SeuratObject) {
  clusters <- unique(SeuratObject$seurat_clusters)
  cluster_info <- data.frame()
  
  for (i in clusters) {
    cluster_data <- subset(SeuratObject, seurat_clusters == i)
    cluster_info <- rbind(cluster_info, data.frame(
      cluster = paste0("Cluster_", i),
      num_of_cells = ncol(cluster_data),
      percentage = signif(ncol(cluster_data) / ncol(SeuratObject) * 100, 3),
      mean_nFeature_RNA = mean(cluster_data@meta.data$nFeature_RNA),
      mean_nCount_RNA = mean(cluster_data@meta.data$nCount_RNA),
      mean_percent.mt = mean(cluster_data@meta.data$percent.mt)
    ))
  }
  
  cluster_info <- rbind(cluster_info, data.frame(
    cluster = "Total",
    num_of_cells = ncol(SeuratObject),
    percentage = 100,
    mean_nFeature_RNA = NA,
    mean_nCount_RNA = NA,
    mean_percent.mt = NA
  ))
  
  table_grob <- tableGrob(cluster_info)
  grid.newpage()
  grid.draw(table_grob)
}

# Load and initialize datasets
datasets <- list(
  "2020_024" = "C:/Users/domin/klaudia_r/Data/2020_024/filtered_feature_bc_matrix",
  "2020_025" = "C:/Users/domin/klaudia_r/Data/2020_025/filtered_feature_bc_matrix",
  "2024_024" = "C:/Users/domin/klaudia_r/Data/2024_024/filtered_feature_bc_matrix",
  "2024_025" = "C:/Users/domin/klaudia_r/Data/2024_025/filtered_feature_bc_matrix",
  "2024_026" = "C:/Users/domin/klaudia_r/Data/2024_026/filtered_feature_bc_matrix",
  "2024_027" = "C:/Users/domin/klaudia_r/Data/2024_027/filtered_feature_bc_matrix"
)

# Create Seurat objects for each dataset
for (name in names(datasets)) {
  assign(paste0("LPS_", name), {
    data <- Read10X(datasets[[name]])
    obj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 500)
    obj <- AddMetaData(obj, metadata = name, col.name = 'Sample_ID')
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
    obj
  })
}

# Check cell counts before and after Seurat object creation
for (name in names(datasets)) {
  raw_data <- get(paste0("LPS_", name, ".data"))
  obj <- get(paste0("LPS_", name))
  cat(paste0("Number of cells in raw data ", name, ": "), ncol(raw_data), "\n")
  cat(paste0("Number of cells in SeuratObject ", name, ": "), ncol(obj), "\n")
}

# Data filtering for LPS_20_024 - The same for each sample
LPS_20_024 <- get("LPS_20_024")
sds_nFeature_RNA <- get_median_sds_as_v(LPS_20_024$nFeature_RNA, TRUE)
sds_nCount_RNA <- get_median_sds_as_v(LPS_20_024$nCount_RNA, FALSE)
sds_percent.mt <- get_median_sds_as_v(LPS_20_024$percent.mt, FALSE)

# Violin plots for initial quality check
colours <- c("red", "orange", "yellow", "green", "yellow", "orange", "red")
colors_two_sd <- c("red", "orange", "yellow", "greenyellow", "green", "greenyellow", "yellow", "orange", "red")
p1 <- VlnPlot(LPS_20_024, features = "nFeature_RNA", group.by = "Sample_ID") + 
  geom_hline(yintercept = sds_nFeature_RNA, colour = colors_two_sd) +
  theme(legend.position = "none")
p2 <- VlnPlot(LPS_20_024, features = "nCount_RNA", group.by = "Sample_ID") + 
  geom_hline(yintercept = sds_nCount_RNA, colour = colours) +
  theme(legend.position = "none")
p3 <- VlnPlot(LPS_20_024, features = "percent.mt", group.by = "Sample_ID") + 
  geom_hline(yintercept = sds_percent.mt, colour = colours) +
  theme(legend.position = "none")
p1 + p2 + p3 + plot_layout(ncol = 3)

# Filter cells based on threshold values
LPS_20_024 <- subset(LPS_20_024, subset = nFeature_RNA < sds_nFeature_RNA[7] & percent.mt < sds_percent.mt[5] & nCount_RNA < sds_nCount_RNA[5])
LPS_20_024 <- subset(LPS_20_024, subset = nFeature_RNA >= 500 & nCount_RNA >= 1000)

# Normalization and clustering
LPS_20_024 <- SCTransform(LPS_20_024, vars.to.regress = "percent.mt", verbose = FALSE)
LPS_20_024 <- RunPCA(LPS_20_024, verbose = FALSE)

# Determine optimal number of PCs
pct <- LPS_20_024[["pca"]]@stdev / sum(LPS_20_024[["pca"]]@stdev) * 100
co1 <- which(cumsum(pct) > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
ElbowPlot(LPS_20_024, ndims = max(30, co1, co2)) + geom_vline(xintercept = c(co1, co2), colour = c("red", "orange"))

# UMAP and clustering
LPS_20_024 <- RunUMAP(LPS_20_024, dims = 1:30, verbose = FALSE)
DimPlot(LPS_20_024, group.by = "Sample_ID") + labs(title = "30 PCs") + theme(legend.position = "none")

# Find clusters at various resolutions
LPS_20_024 <- FindNeighbors(LPS_20_024, dims = 1:30, verbose = FALSE)
for (res in c(0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4)) {
  LPS_20_024 <- FindClusters(LPS_20_024, resolution = res, algorithm = 1, verbose = FALSE)
}
clustree(LPS_20_024@meta.data, prefix = "SCT_snn_res.")

# Set clusters and plot UMAP
LPS_20_024$seurat_clusters <- LPS_20_024$SCT_snn_res.1
Idents(LPS_20_024) <- LPS_20_024$SCT_snn_res.1
DimPlot(LPS_20_024, label = TRUE)

# Generate cluster information table
tabulate_seurat_clusters_info(LPS_20_024)

# Save Seurat object
saveRDS(LPS_20_024, "C:/Users/domin/klaudia_r/Data_BetterQuality/2020_024/LPS_20_024.rds")
