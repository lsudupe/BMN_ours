## SCRIPT: Spatial object clustering to check heterogeneity and 
## analize strange cluster BM project

## 23.12.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2) 
library(tidyverse)
library(base)


## Spatial object
M1_fem_1C <- readRDS("./objects/card/M1_fem_1C_subgroup.rds")
M1_tib_1A <- readRDS("./objects/card/M1_tib_1A_subgroup.rds")
M3_fem_1C <- readRDS("./objects/card/M3_fem_1C_subgroup.rds")
M3_tib_2A <- readRDS("./objects/card/M3_tib_2A_subgroup.rds")

##Merge them
combined <- merge(M1_fem_1C, y = c( M3_fem_1C), 
                  add.cell.ids = c("M1_fem_1C", "M3_fem_1C"), project = "BM")

##Merge them
#combined <- merge(M1_tib_1A, y = c( M3_tib_2A), 
#                  add.cell.ids = c("M1_tib_1A", "M3_tib_2A"), project = "BM")

##Merge them
#combined <- merge(M1_tib_1A, y = c( M3_tib_2A, M1_fem_1C, M3_fem_1C), 
#                  add.cell.ids = c("M1_tib_1A", "M3_tib_2A", "M1_fem_1C", "M3_fem_1C"), project = "BM")

x <- combined
x.image <- x@images
x@images[["M1_fem_1C"]]<- NULL
x@images[["M3_fem_1C"]]<- NULL

#x@images[["M1_tib_1A"]]<- NULL
#x@images[["M3_tib_2A"]]<- NULL

######seurat#####
list <- SplitObject(x, split.by = "type")
list <- lapply(X = list, FUN = SCTransform, assay="Spatial")
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
x <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

x <- RunPCA(x,npcs = 20, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution=0.8) 

x@images <- x.image

samples <- c(x)
names(samples) <- c("femur")
#names(samples) <- c("tibia")
#names(samples) <- c("all")

for (i in 1:length(samples)){
  a <- samples[[i]]
  #umap separate in cell type
  #pdf(file.path("./results/heterogeneity/femur/",filename = paste0("umap_clusters_",names(samples[i]),".pdf")))
  #pdf(file.path("./results/heterogeneity/tibia/",filename = paste0("umap_clusters_",names(samples[i]),".pdf")))
  pdf(file.path("./results/heterogeneity/all/",filename = paste0("umap_clusters_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("seurat_clusters"), label = T) + ggtitle("cell type"))
  dev.off()
  #umap separate in samples
  #pdf(file.path("./results/heterogeneity/femur/",filename = paste0("umap_samples_",names(samples[i]),".pdf")))
  #pdf(file.path("./results/heterogeneity/tibia/",filename = paste0("umap_samples_",names(samples[i]),".pdf")))
  pdf(file.path("./results/heterogeneity/all/",filename = paste0("umap_samples_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("orig.ident"), label = T) + ggtitle("sample"))
  dev.off()
  #umap separate in samples
  #pdf(file.path("./results/clusters/integration/",filename = paste0("umap_area_",names(samples[i]),".pdf")))
  #print(DimPlot(a, group.by = c("area"),repel=TRUE, label = T) + ggtitle("sample"))
  #dev.off()
  #umap separate in disease
  #pdf(file.path("./results/heterogeneity/femur/",filename = paste0("umap_type_",names(samples[i]),".pdf")))
  #pdf(file.path("./results/heterogeneity/tibia/",filename = paste0("umap_type_",names(samples[i]),".pdf")))
  pdf(file.path("./results/heterogeneity/all/",filename = paste0("umap_type_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("type"), label = T) + ggtitle("type"))
  dev.off()
  a <- SetIdent(a, value = a@meta.data[["seurat_clusters"]])
  #spatial umap 
  #pdf(file.path("./results/heterogeneity/femur/",filename = paste0("cluster_spatial_",names(samples[i]),".pdf")))
  #pdf(file.path("./results/heterogeneity/tibia/",filename = paste0("cluster_spatial_",names(samples[i]),".pdf")))
  pdf(file.path("./results/heterogeneity/all/",filename = paste0("cluster_spatial_",names(samples[i]),".pdf")))
  print(SpatialDimPlot(a, combine = FALSE,label.size = 1.5, label = T, crop = TRUE, pt.size.factor = 10))
  dev.off()
  #a <- SetIdent(a, value = a@meta.data[["area"]])
  #spatial umap 
  #pdf(file.path("./results/clusters/integration/",filename = paste0("area_spatial_",names(samples[i]),".pdf")))
  #print(SpatialDimPlot(a, combine = FALSE,label.size = 1.5, label = T))
  #dev.off()
}


##################################CELL TYPE PROPORTIONS HETEROGENEITY
single_cell_bonemarrow <- readRDS("./objects/heterogeneity/single_cell_bonemarrow.rds")







##################################AMAIA PART HETEROGENEITY

#ERYTHROBLAST. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7689499/
CCL3 #up
KLF1 #down
GATA1 #down
# Osteoblast. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4580250/
DKK1 # up
MIP-1
CCR1 #which is one pathway contributing to their decreased osteogenic capacity
CXCL12



DefaultAssay(x) <- "SCT"
DefaultAssay(x) <- "integrated"
SpatialDimPlot(x, features = c("active.ident"))
FeaturePlot(x, features = "Ccr1")
DimPlot(x, group.by = c("seurat_clusters"), label = T)
SpatialDimPlot(x, label = T, crop = TRUE, pt.size.factor = 10)
SpatialFeaturePlot(x,  features = c("Cxcl12"), pt.size = 10)
SpatialFeaturePlot(x, features = c("Klf1"), pt.size = 10)
SpatialFeaturePlot(x, features = c("Ccl3"), pt.size = 10)

markers_seurat_area <- Seurat::FindAllMarkers(object = x, 
                                              assay = "integrated",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)

markers_seurat <- subset(markers_seurat_area, p_val_adj < 0.05 & 0.5 < avg_log2FC)

#Top5
top20_area <- markers_seurat %>%
  group_by(cluster) %>%
  top_n(n = 20,
        wt = avg_log2FC)

DoHeatmap(x, features = top20_area$gene) + 
  theme(text = element_text(size = 4.5))

cell_types <- as.list(as.vector(unique(top20_area$cluster)))
names(cell_types) <- as.vector(unique(top20_area$cluster))

for (i in 1:length(cell_types)){
  a <- cell_types[[i]]
  top20_area$cluster <- as.character(top20_area$cluster)
  b <- top20_area[top20_area$cluster %in% c(a), ]
  b <- as.vector(b$gene)
  cell_types[[i]] <- b
}

cluster_5 <- as.vector(cell_types[6])

DefaultAssay(x) <- "SCT"
DefaultAssay(x) <- "integrated"
DefaultAssay(x) <- "Spatial"
SpatialFeaturePlot(x, features = c("Bst2"), pt.size = 10)


p1 <- SpatialFeaturePlot(spatial, features = "Fcgr1", pt.size.factor = 3, alpha= 0.6, combine = FALSE)
fix.p1 <- scale_fill_gradientn(limits = c(0,15),
                               breaks=c("Min","Max"),
                               labels=c("Min","Max"),
                               colours=topo.colors(7))
p2 <- lapply(p1, function (x) x + fix.p1)

pdf("./f.pdf")
CombinePlots(p2)
dev.off()


p1 <- SpatialFeaturePlot(spatial, features = "Fcgr1", pt.size.factor = 3, alpha = 0.6, combine = FALSE)
fix.p1 <- scale_fill_continuous(type = "viridis")
p2 <- lapply(p1, function (x) x + fix.p1)

pdf(paste(i,".pdf",sep=""))
print(CombinePlots(p2))
dev.off()

##################################AMAIA PART HETEROGENEITY
