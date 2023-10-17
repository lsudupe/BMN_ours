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
  
  ##ratio
  enrich7 <- a@meta.data[["signature_1_cluster7"]]
  enrich4 <- a@meta.data[["signature_1_cluster4"]]
  enrich1 <- a@meta.data[["signature_1_cluster1"]]
  
  #ratio 1 vs 7
  ratio_cluster1vs7 <- log10(enrich1/enrich7)
  is.na(ratio_cluster1vs7) <-sapply(ratio_cluster1vs7, is.infinite)
  ratio_cluster1vs7[is.na(ratio_cluster1vs7)] = 0
  a$ratio_cluster1vs7 <- ratio_cluster1vs7
  
  #plot
  bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
  re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
  
  l <- c(min(a@meta.data[["ratio_cluster1vs7"]]), max(a@meta.data[["ratio_cluster1vs7"]]))
  p1 <- SpatialFeaturePlot(a, features = c("ratio_cluster1vs7"), combine = FALSE, ncol = 2)
  fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re),
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 na.value = "grey98",
                                 limits = l)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0("./results/human/cluster_enrich/ratio/", b, "_", "ratio_cluster1vs7.pdf", sep = ""))
  print(CombinePlots(p2))
  dev.off()
  
  #ratio 4 vs 7
  ratio_cluster4vs7 <- log10(enrich4/enrich7)
  is.na(ratio_cluster4vs7) <-sapply(ratio_cluster4vs7, is.infinite)
  ratio_cluster4vs7[is.na(ratio_cluster4vs7)] = 0
  a$ratio_cluster4vs7 <- ratio_cluster4vs7
  
  #plot
  bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
  re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
  
  l <- c(min(a@meta.data[["ratio_cluster4vs7"]]), max(a@meta.data[["ratio_cluster4vs7"]]))
  p1 <- SpatialFeaturePlot(a, features = c("ratio_cluster4vs7"), combine = FALSE, ncol = 2)
  fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re),
                                 breaks=b,
                                 labels=c("Min","Max"),
                                 na.value = "grey98",
                                 limits = l)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(paste0("./results/human/cluster_enrich/ratio/", b, "_", "ratio_cluster4vs7.pdf", sep = ""))
  print(CombinePlots(p2))
  dev.off()
  
  
  post[[length(post) + 1]] <- a
}


