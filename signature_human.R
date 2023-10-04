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


Idents(se) <- "clustering"
all_markers <- FindAllMarkers(se, min.pct = 0.25)

top_markers_per_cluster <- all_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 25, wt = -log10(p_val_adj))

write.csv(top_markers_per_cluster, file = "./25.csv", row.names = FALSE)

split_data <- split(top_markers_per_cluster, top_markers_per_cluster$cluster)

# extract genes per cluster
for(i in 1:length(split_data)) {
  cluster_name <- unique(split_data[[i]]$cluster)[1]
  gene_list <- toupper(split_data[[i]]$gene)  # Convert gene names to uppercase
  assign(paste0("cluster", cluster_name), gene_list)
}

# split human
B08041 <- human[["BM_human_AP-B08041_"]]
B08805 <- human[["BM_human_AP-B08805"]]
B10395 <- human[["BM_B10395"]]

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
    
    file_name <- paste0("./results/human/cluster_enrich/new/100/", b, "_", cluster_name, "_pc.pdf", sep = "")
    
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

  post[[length(post) + 1]] <- a
}


