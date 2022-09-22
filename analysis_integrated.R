## SCRIPT: Integration spatial data BMN project

## 20.09.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2) 
library(harmony)
library('magrittr')
library(tidyverse)
library(base)
library(scclusteval)

#Data--------------------------------------
combined  <- readRDS("./objects/sp/combined_filtered.rds")
x <- combined

x.image <- x@images
x@images[["M1_tib_1A"]]<- NULL
x@images[["M1_fem_2B"]]<- NULL
x@images[["M3_fem_1C"]]<- NULL
x@images[["M3_tib_1A"]]<- NULL

###separate the data
seurat_resolution <- x
harmony_resolution <- x

######seurat#####
list <- SplitObject(seurat_resolution, split.by = "type")
list <- lapply(X = list, FUN = SCTransform, assay="Spatial")
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
seurat_resolution <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

seurat_resolution <- RunPCA(seurat_resolution,npcs = 15, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:15, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:15) %>%
  FindClusters(resolution=0.6) 

seurat_resolution@images <- x.image
saveRDS(seurat_resolution, "./objects/sp/integrated/integrated.seurat_type.rds")

######harmony#####
harmony_resolution <- harmony_resolution %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(npcs = 15) %>%
  RunHarmony(assay.use="Spatial",reduction = "pca", dims = 1:15, group.by.vars = "type") %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution=0.6) %>%
  RunUMAP(reduction = "harmony", dims = 1:15, n.epochs = 1e3) 

harmony_resolution@images <- x.image
saveRDS(harmony_resolution, "./objects/sp/integrated/integrated.harmony_type.rds")

integrated_seurat <- readRDS("./objects/sp/integrated/integrated.seurat_type.rds")
integrated_harmony <- readRDS("./objects/sp/integrated/integrated.harmony_type.rds")

######Jackard index
pdf(file.path("./results/jackard",filename = "jackard_type.pdf"))                                                         
PairWiseJaccardSetsHeatmap(integrated_seurat@active.ident,
                           integrated_harmony@active.ident)
dev.off()

samples <- c(integrated_seurat, integrated_harmony)
names(samples) <- c("integrated_seurat", "integrated_harmony")

for (i in 1:length(samples)){
  a <- samples[[i]]
  #umap separate in cell type
  pdf(file.path("./results/clusters/integration/",filename = paste0("type_bycelltype",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("seurat_clusters"), label = T) + ggtitle("cell type"))
  dev.off()
  #umap separate in samples
  pdf(file.path("./results/clusters/integration/",filename = paste0("type_bysamples",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("orig.ident"), label = T) + ggtitle("sample"))
  dev.off()
  #umap separate in samples
  pdf(file.path("./results/clusters/integration/",filename = paste0("type_bytype",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("type"), label = T) + ggtitle("sample"))
  dev.off()
  #spatial umap 
  pdf(file.path("./results/clusters/integration/",filename = paste0("type_spatial",names(samples[i]),".pdf")))
  print(SpatialDimPlot(a, image.alpha = 1, alpha = 1, combine = FALSE, label = T))
  dev.off()
  
}
