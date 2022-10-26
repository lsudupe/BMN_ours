## SCRIPT: Check cell type heterogeneity in my areas BM project

## 26.10.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)

#Data---------------------------------
spatial <- readRDS("./objects/sp/second/combined_filtered.rds")

######seurat#####
x <- spatial

x.image <- x@images
x@images[["M1_tib_1A"]]<- NULL
x@images[["M1_fem_1C"]]<- NULL
x@images[["M3_tib_2A"]]<- NULL
x@images[["M3_fem_1C"]]<- NULL

###separate the data
seurat_resolution <- x

######seurat#####
list <- SplitObject(seurat_resolution, split.by = "orig.ident")
list <- lapply(X = list, FUN = SCTransform, assay="Spatial")
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
seurat_resolution <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

seurat_resolution@images <- x.image
seurat_resolution <- RunPCA(seurat_resolution,npcs = 15, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:15, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:15) %>%
  FindClusters(resolution=0.7) 

seurat_resolution <- SetIdent(seurat_resolution, value = seurat_resolution@meta.data[["area"]])
#spatial areas 
pdf(file.path("./results/area_manual/",filename = "area_spatial_.pdf"))
print(SpatialDimPlot(seurat_resolution, combine = FALSE,label = T))
dev.off()


#########################DE#########################3
seurat_resolution <- SetIdent(seurat_resolution, value = seurat_resolution@meta.data[["type"]])
markers <- Seurat::FindAllMarkers(object = seurat_resolution, 
                                      assay = "integrated",
                                      verbose = TRUE)

markers_conserved <- FindConservedMarkers(object = seurat_resolution, 
                                  assay = "integrated",
                                  ident.1= "MM",
                                  ident.2= "control",
                                  verbose = TRUE)

markers_area <- subset(markers, p_val_adj < 0.05 & 0.5 < avg_log2FC)
top20_area <- markers %>%
  group_by(cluster) %>%
  top_n(n = 200,
        wt = avg_log2FC)


pdf(file.path("./results/area_manual/",filename = "heatmap.pdf"))
print(DoHeatmap(seurat_resolution, assay= "integrated",features = top20_area$gene))
dev.off()


