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

## analysis
x <- combined
images <- x@images

## samples
samples <- unique(combined@meta.data[["orig.ident"]])
for (i in samples){
  x@images[[i]] <- NULL
}

x <- SCTransform(x, assay = "Spatial", verbose = TRUE, method = "poisson")
x <- RunPCA(x, assay = "SCT", verbose = FALSE)
x <- FindNeighbors(x, reduction = "pca", dims = 1:30)
x <- FindClusters(x, verbose = FALSE)
x <- RunUMAP(x, reduction = "pca", dims = 1:30)

Seurat::Idents(object = x) <- combined@meta.data[["type"]]

###separate the data
seurat_resolution <- x
harmony_resolution <- x

######seurat#####
#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Dimensionality_reduction_and_clustering
# need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 2000 * 1024^2)  # set allowed size to 2K MiB

list <- SplitObject(seurat_resolution, split.by = "orig.ident")
list <- lapply(X = list, FUN = SCTransform, assay="Spatial", method = "poisson")
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000, verbose = FALSE)
list <- PrepSCTIntegration(object.list = list, anchor.features = features, verbose = FALSE)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features, verbose = FALSE)
integrated_seurat <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

seurat_resolution <- RunPCA(integrated_seurat, verbose = FALSE) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution=0.7) %>%
  RunUMAP(dims = 1:20) 

seurat_resolution@images <- images
saveRDS(seurat_resolution, "./objects/sp/integrated/integrated.seurat.rds")


######harmony#####
harmony_resolution <- harmony_resolution %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(npcs = 20) %>%
  RunHarmony(assay.use="Spatial",reduction = "pca", dims = 1:20, group.by.vars = "orig.ident") %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution=0.7) %>%
  RunUMAP(reduction = "harmony", dims = 1:20, n.epochs = 1e3) 

harmony_resolution@images <- images
saveRDS(harmony_resolution, "./objects/sp/integrated/integrated.harmony.rds")

integrated_seurat <- readRDS("./objects/sp/integrated/integrated.seurat.rds")
integrated_harmony <- readRDS("./objects/sp/integrated/integrated.harmony.rds")

######Jackard index harmony clusters vs seurat
pdf(file.path("./results/jackard/",filename = "jackard_type_harmonyvsseurat.pdf"))                                                         
PairWiseJaccardSetsHeatmap(integrated_seurat@active.ident,
                           integrated_harmony@active.ident,
                           best_match = TRUE)
dev.off()

## spatial plots

samples <- c(integrated_seurat, integrated_harmony)
names(samples) <- c("integrated_seurat", "integrated_harmony")


for (i in 1:length(samples)){
  a <- samples[[i]]
  #umap separate in cell type
  pdf(file.path("./results/clusters/integration/orig.ident",filename = paste0("umap_clusters_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("seurat_clusters"), label = T) + ggtitle("cell type"))
  dev.off()
  #umap separate in samples
  pdf(file.path("./results/clusters/integration/orig.ident",filename = paste0("umap_samples_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("orig.ident"), label = T) + ggtitle("sample"))
  dev.off()
  #umap separate in samples
  #pdf(file.path("./results/clusters/integration/",filename = paste0("umap_area_",names(samples[i]),".pdf")))
  #print(DimPlot(a, group.by = c("area"),repel=TRUE, label = T) + ggtitle("sample"))
  #dev.off()
  #umap separate in samples
  pdf(file.path("./results/clusters/integration/orig.ident",filename = paste0("umap_type_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("type"), label = T) + ggtitle("sample"))
  dev.off()
  a <- SetIdent(a, value = a@meta.data[["seurat_clusters"]])
  #spatial umap 
  pdf(file.path("./results/clusters/integration/orig.ident",filename = paste0("cluster_spatial_",names(samples[i]),".pdf")))
  print(SpatialDimPlot(a, combine = FALSE,label.size = 1.5, crop=TRUE,label = T, pt.size.factor = 7))
  dev.off()
  #a <- SetIdent(a, value = a@meta.data[["area"]])
  #spatial umap 
  #pdf(file.path("./results/clusters/integration/",filename = paste0("area_spatial_",names(samples[i]),".pdf")))
  #print(SpatialDimPlot(a, combine = FALSE,label.size = 1.5, label = T))
  #dev.off()
}

## check samples in object
sample <- c(unique(integrated_seurat$orig.ident))

## take out image
images <- integrated_seurat@images

## set all the images null
for (i in sample){
  integrated_seurat@images[[i]] <- NULL
}

## separate samples in list
list <- SplitObject(integrated_seurat, split.by = "orig.ident")

## add to each element in list its image
samples <- c()
for (i in sample){
  a <- list[[i]]
  a@images[[i]] <- images[[i]]
  samples[[length(samples) + 1]] <- a
}
names(samples) <- sample

## separate list
#list2env(list,envir=.GlobalEnv)

for (i in 1:length(samples)){
  a <- samples[[i]]
  pdf(file.path("./results/clusters/individual/",filename = paste0("cluster_spatial_",names(samples[i]),".pdf")))
  print(SpatialDimPlot(a ,label = T,crop = TRUE, pt.size.factor = 7))
  dev.off()
}

x <- integrated_seurat

###Markers
x <- integrated_seurat
x <- SetIdent(x, value = x@meta.data[["area"]])
markers_seurat_area <- Seurat::FindAllMarkers(object = x, 
                                    assay = "Spatial",
                                    slot = "data",
                                    verbose = TRUE, 
                                    only.pos = TRUE)

x <- SetIdent(x, value = area@meta.data[["seurat_clusters"]])
markers_seurat_clusters <- Seurat::FindAllMarkers(object = x, 
                                         assay = "Spatial",
                                         slot = "data",
                                         verbose = TRUE, 
                                         only.pos = TRUE)

#saveRDS(markers_x, "./results/marker_genes/M1_tib_1A.rds")

#Filter
markers_area <- subset(markers_seurat_area, p_val_adj < 0.05 & 0.5 < avg_log2FC)
markers_seurat <- subset(markers_seurat_clusters, p_val_adj < 0.05 & 0.5 < avg_log2FC)

#Top5
top20_area <- markers_area %>%
  group_by(cluster) %>%
  top_n(n = 20,
        wt = avg_log2FC)

top20_seurat <- markers_seurat %>%
  group_by(cluster) %>%
  top_n(n = 20,
        wt = avg_log2FC)

write.csv(top20_area, "./results/DE/emma/top20_genes_area.csv", row.names=FALSE)


cell_types <- as.list(as.vector(unique(top20_area$cluster)))
names(cell_types) <- as.vector(unique(top20_area$cluster))

for (i in 1:length(cell_types)){
  a <- cell_types[[i]]
  top20_area$cluster <- as.character(top20_area$cluster)
  b <- top20_area[top20_area$cluster %in% c(a), ]
  b <- as.vector(b$gene)
  cell_types[[i]] <- b
}

integrated_seurat <- readRDS("./objects/sp/integrated/integrated.seurat_type.rds")


####AddmoduleScore with Mahtab annotation
DefaultAssay(integrated_seurat) <- "SCT"
for (i in 1:length(cell_types)){
  a <- cell_types[[i]]
  c <- paste0(names(cell_types[i]), "1")
  integrated_seurat <- AddModuleScore(integrated_seurat, features = cell_types[i], name = names(cell_types[i]))
  #spatial 
  pdf(file.path("./results/marker_genes/areas_Mahtab/seurat/",filename = paste0("spatial_seurat_",c,".pdf")))
  print(SpatialFeaturePlot(integrated_seurat,features=c,combine = FALSE))
  dev.off()
}

DoHeatmap(integrated_seurat,assay = "SCT", group.by = integrated_seurat@meta.data[["type"]],
          features = cell_types[[1]])
DoHeatmap(integrated_seurat,assay = "SCT", features = top20_seurat$gene)
