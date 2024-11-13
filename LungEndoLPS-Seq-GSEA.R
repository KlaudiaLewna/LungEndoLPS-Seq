# Load necessary libraries
library(dplyr)
library(fgsea)
library(clusterProfiler)
library(msigdbr)
library(DOSE)  # For enrichment analysis
library(enrichplot)  # For visualization
library(data.table)
library(ggplot2)
library(forcats)
library(purrr)

# Set working directory and seed
setwd("C:\\Users\\domin\\klaudia_r\\DE_analysis_scrnaseq\\Results_Harmony")
set.seed(43)

# Define the cluster names and corresponding gene file names
clusters <- c("aCap", "gCap", "Vein", "Artery", "LEC", "gCap-H2", "cap-Fabp4+", "cap-Lysmd4+", "cap-Tubb2b")
gene_files <- paste0(clusters, "_signif_genes.csv")

# Initialize a list to store ranked gene lists
ranked_gene_lists <- list()

# Loop through clusters to generate ranked gene lists
for (i in seq_along(clusters)) {
  # Read the significant gene data for the current cluster
  significant_genes <- read.table(gene_files[i], header = TRUE, sep = ",")
  
  # Filter out rows where log2FoldChange is NA
  gseaDat <- filter(significant_genes, !is.na(log2FoldChange))
  
  # Create a rank list based on log2FoldChange and store it
  ranks <- gseaDat$log2FoldChange
  names(ranks) <- gseaDat$gene
  ranked_gene_lists[[clusters[i]]] <- ranks
}

# Load hallmark pathways (or other pathway sets)
hallmark_pathways <- gmtPathways("mh.all.v2024.1.Mm.symbols.gmt")

# Initialize list for storing GSEA results
gsea_results_list <- list()

# Perform GSEA for each cluster
for (i in seq_along(clusters)) {
  cluster_name <- clusters[i]
  ranks <- ranked_gene_lists[[cluster_name]]
  
  # Run GSEA
  fgseaRes <- fgsea(
    pathways = hallmark_pathways, 
    stats = ranks, 
    minSize = 15, 
    maxSize = 500, 
    nperm = 1000
  )
  
  fgseaRes <- fgseaRes %>%
    mutate(Cluster = cluster_name)  # Add cluster name
  
  # Store GSEA results for each cluster
  gsea_results_list[[cluster_name]] <- fgseaRes
}

# Combine all GSEA results into one data frame
combined_gsea_results <- bind_rows(gsea_results_list)

# KEGG Pathways: Load KEGG pathways for mouse from MSigDB
msigdb_mouse <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
kegg_pathways_list <- split(msigdb_mouse$gene_symbol, msigdb_mouse$gs_name)

# Perform GSEA for KEGG pathways
gsea_results_list_kegg <- list()
for (cluster_name in names(ranked_gene_lists)) {
  cluster_ranks <- ranked_gene_lists[[cluster_name]]
  
  # Run GSEA for KEGG pathways
  gsea_res <- fgsea(
    pathways = kegg_pathways_list,
    stats = cluster_ranks,
    minSize = 15,
    maxSize = 500,
    nperm = 1000
  )
  
  gsea_res$Cluster <- cluster_name
  gsea_results_list_kegg[[cluster_name]] <- gsea_res
}

# Combine all KEGG GSEA results into one data frame
combined_gsea_results_kegg <- bind_rows(gsea_results_list_kegg)

# Filter GSEA results for selected pathways
dot_df_filtered <- combined_gsea_results_kegg %>%
  filter(pathway %in% selected_pathways)

# Top pathways based on NES
dot_df_filtered <- combined_gsea_results_kegg %>%
  arrange(desc(abs(NES))) %>%
  top_n(210, wt = abs(NES))

# Plot the dot plot with custom color scale for padj
ggplot(dot_df_filtered, aes(x = Cluster, y = fct_reorder(pathway, NES))) + 
  geom_point(aes(size = NES, color = padj)) +
  scale_size_continuous(range = c(0, 5)) +
  scale_colour_gradientn(colours = c("green", "purple"), 
                         values = scales::rescale(c(0, 0.05, 0.2)), 
                         limits = c(0, 0.2)) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    panel.spacing = unit(1, "lines") 
  ) +
  ylab(NULL) +
  xlab("Cluster")
