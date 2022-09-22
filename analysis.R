## SCRIPT: Cluster analysis of the spatial data BMN project

## 15.09.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2) 
library('magrittr')
library(tidyverse)
library(base)

#Data--------------------------------------
M1_tib_1A  <- readRDS("./objects/sp/M1_tib_1A_filtered.rds")
x <- M1_tib_1A

x <- SCTransform(x, assay="Spatial")
x <- RunPCA(x, verbose = FALSE)
x <- RunUMAP(x, reduction = "pca", dims = 1:20, verbose = FALSE)
x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
x <- FindClusters(x, resolution=0.8)

x@meta.data[["orig.ident"]] <- "M1_tib_1A"
saveRDS(x, "./objects/sp/M1_tib_1A_normalized.rds")
write.csv(x@meta.data,"./objects/sp/M1_tib_1A_normalized_meta.csv", row.names = FALSE)

#umap separate in cell type
pdf(file.path("./results/clusters",filename = "bycelltype_M1_tib_1A.pdf"))
print(DimPlot(x, group.by = c("seurat_clusters"), label = T) + ggtitle("cell type"))
dev.off()

#spatial umap 
pdf(file.path("./results/clusters",filename = "bycelltype_spatial_M1_tib_1A.pdf"))
print(SpatialDimPlot(x, label = T))
dev.off()

x <- readRDS("./objects/sp/M1_tib_1A_normalized.rds")
###Markers

markers_x <- Seurat::FindAllMarkers(object = x, 
                                          assay = "SCT",
                                          slot = "data",
                                          verbose = TRUE, 
                                          only.pos = TRUE)

saveRDS(markers_x, "./results/marker_genes/M1_tib_1A.rds")

#Filter
markers_x <- subset(markers_x, p_val_adj < 0.05 & 0.5 < avg_log2FC)

#Top5
top3_x <- markers_x %>%
  group_by(cluster) %>%
  top_n(n = 3,
        wt = avg_log2FC)

####Markers per cluster
top3_x$cluster <- as.character(top3_x$cluster)
cluster_0 <- top3_x[top3_x$cluster %in% c("0"), ]
cluster_0 <- as.vector(cluster_0$gene)

cluster_1 <- top3_x[top3_x$cluster %in% c("1"), ]
cluster_1 <- as.vector(cluster_1$gene)

cluster_2 <- top3_x[top3_x$cluster %in% c("2"), ]
cluster_2 <- as.vector(cluster_2$gene)

cluster_3 <- top3_x[top3_x$cluster %in% c("3"), ]
cluster_3 <- as.vector(cluster_3$gene)

cluster_4 <- top3_x[top3_x$cluster %in% c("4"), ]
cluster_4 <- as.vector(cluster_4$gene)

cluster_5 <- top3_x[top3_x$cluster %in% c("5"), ]
cluster_5 <- as.vector(cluster_5$gene)

cluster_6 <- top3_x[top3_x$cluster %in% c("6"), ]
cluster_6 <- as.vector(cluster_6$gene)


list_top <- list(cluster_0, cluster_1, cluster_2, cluster_3, cluster_4, cluster_5, cluster_6)
names(list_top) <- c("cluster_0", "cluster_1", "cluster_2","cluster_3", "cluster_4", "cluster_5", "cluster_6")


for (i in 1:length(list_top)){
  a <- list_top[[i]]
  p1 <- SpatialFeaturePlot(x, features = a, ncol = 3 )
  pdf(file.path("./results/marker_genes/",filename = paste0(i,".spatial_genes.pdf",sep="")))
  print(p1)
  dev.off()
  p1 <- FeaturePlot(x, features = a, pt.size = 0.2, ncol = 1 )
  pdf(file.path("./results/marker_genes/",filename = paste0(i,".feature_genes.pdf",sep="")))
  print(p1)
  dev.off()
  
}





