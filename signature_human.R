## SCRIPT: Extract signature per cluster and check in human BM project

## 03.10.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(RColorBrewer)
library(tidyverse)
library(UCell)

#Data---------------------------------
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")
human <- readRDS("./objects/sp/human/human_combined.rds")

# split human
B08041 <- human[["BM_human_AP-B08041_"]]
B08805 <- human[["BM_human_AP-B08805"]]
B10395 <- human[["BM_B10395"]]

Idents(se) <- "clustering"
all_markers <- FindAllMarkers(se, min.pct = 0.25)

top_markers_per_cluster <- all_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = -log10(p_val_adj))

split_data <- split(top_markers_per_cluster, top_markers_per_cluster$cluster)

# Assume you have gene lists for each human sample
genes_B08041 <- rownames(B08041)  # Extract gene names
genes_B08805 <- rownames(B08805)
genes_B10395 <- rownames(B10395)

# Find common genes across all three samples
common_genes <- Reduce(intersect, list(genes_B08041, genes_B08805, genes_B10395))

# Function to get top genes ensuring they are in the common gene list
get_top_genes <- function(cluster_data, common_genes, n = 25) {
  # Convert gene names in the cluster data to uppercase
  cluster_genes_upper <- toupper(cluster_data$gene)
  # Ensure that common_genes are also in uppercase
  common_genes_upper <- toupper(common_genes)
  # Check for common genes in uppercase
  common_cluster_genes <- cluster_genes_upper[cluster_genes_upper %in% common_genes_upper]
  if(length(common_cluster_genes) >= n) {
    return(head(common_cluster_genes, n))
  } else {
    # Select additional genes if fewer than n are common
    additional_genes <- setdiff(cluster_genes_upper, common_cluster_genes)
    total_genes <- c(common_cluster_genes, head(additional_genes, n - length(common_cluster_genes)))
    return(total_genes)
  }
}

# Your existing loop
for(i in 1:length(split_data)) {
  cluster_name <- unique(split_data[[i]]$cluster)[1]
  top_genes <- get_top_genes(split_data[[i]], common_genes)
  assign(paste0("cluster", cluster_name), top_genes)
}

lista <- c(B08041, B08805, B10395)
names(lista) <- c("B08041","B08805", "B10395")

###Enrichment score
color <- brewer.pal(11,"Spectral")
color <- rev(color)

post <- c()
# Define the cluster names
cluster_names <- paste0("cluster", 1:7)

for (i in 1:length(lista)){
  v <- "sss"
  a <- lista[[i]]
  b <- names(lista[i])
  
  for (cluster_name in cluster_names) {
    # Define the variable names dynamically based on the cluster_name
    signature_var <- paste0("signature_1_", cluster_name)
    
    file_name <- paste0("./results/human/cluster_enrich/new/25/", b, "_", cluster_name, "_pc.pdf", sep = "")
    
    # Add UCellScore for the current cluster
    vector <- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(get(cluster_name)))
    a@meta.data[[signature_var]] <- as.vector(vector)
    
    #c <- c(min(a@meta.data[[signature_var]]), max(a@meta.data[[signature_var]]))
    # Create the plot
    p <- FeatureOverlay(a, features = c(signature_var), ncols = 1, pt.size = 1.1, 
                        value.scale = "all" ,cols = color) +
      scale_fill_gradientn(colours = color,
                           breaks = c(0.0, 0.5),
                           labels = c("Min", "Max"),
                           limits = c(0.0, 0.5))
    
    # Save the plot to a PDF file
    pdf(file_name)
    print(p)
    dev.off()
  }
  
  #########
  # Call the function for each comparison
  plot_ratio(a, 1, 7, color, "./results/human/cluster_enrich/25_ratio/")
  plot_ratio(a, 5, 7, color, "./results/human/cluster_enrich/25_ratio/")
  plot_ratio(a, 5, 6, color, "./results/human/cluster_enrich/25_ratio/")
  plot_ratio(a, 6, 7, color, "./results/human/cluster_enrich/25_ratio/")
  
  post[[length(post) + 1]] <- a
}

names(post) <- c("B08041","B08805", "B10395")

