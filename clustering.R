## SCRIPT: Spatial object cluster comparation BM project

## 22.09.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(SingleR)
library(dplyr)
library(ggplot2) 
library('magrittr')
library(tidyverse)
library(base)

#Data--------------------------------------
load(file='./data/single-cell/RNAMagnetDataBundle/NicheMarkers10x.rda')
load(file='./data/single-cell/RNAMagnetDataBundle/NicheData10x.rda')
load(file="./data/single-cell/RNAMagnetDataBundle/NicheDataLCM.rda")
load(file="./data/single-cell/RNAMagnetDataBundle/NicheMetaDataLCM.rda")
single_cell_bonemarrow = UpdateSeuratObject(object = NicheData10x)

single_cell_bonemarrow@meta.data[["cell"]] <- single_cell_bonemarrow@active.ident
#####SUBSET IN BONE, BONE MARROW

pdf(file.path("./results/clusters/single_cell//",filename = "cell_types_singlecell.pdf"))
print(DimPlot(single_cell_bonemarrow, reduction = "tsne", label = TRUE))
dev.off()

integrated_seurat <- readRDS("./objects/sp/integrated/integrated.seurat_type.rds")
integrated_harmony <- readRDS("./objects/sp/integrated/integrated.harmony_type.rds")

# Plot the elbow plot
ElbowPlot(object = integrated_harmony, 
          ndims = 40)
DimPlot(integrated_seurat, reduction = "pca")


DimHeatmap(integrated_seurat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(integrated_seurat, dims = 1:15, cells = 500, balanced = TRUE)

integrated_harmony <- JackStraw(integrated_harmony, num.replicate = 100)
integrated_harmony <- ScoreJackStraw(integrated_harmony, dims = 1:20)
JackStrawPlot(integrated_harmony, dims = 1:15)


#####Single-R
single.sce <- as.SingleCellExperiment(single_cell_bonemarrow)
seurat.sce <- as.SingleCellExperiment(integrated_seurat)
harmony.sce <- as.SingleCellExperiment(integrated_harmony)


###seurat predictions
predictions_seurat <- SingleR(test=seurat.sce, assay.type.test=1, 
                          ref=single.sce, labels=single.sce$cell)

table(predictions_seurat$labels)
a <- as.factor(predictions_seurat@listData[["labels"]])
integrated_seurat@meta.data[["pred"]] <- a
b <- integrated_seurat
b@active.ident <- b@meta.data[["pred"]]
saveRDS(b, "./objects/sp/integrated/integrated.seurat_predictions.rds")

pdf(file.path("./results/singleR/",filename = "cell_types_singlecell_seurat.pdf"))
print(SpatialDimPlot(b, combine = FALSE,label=TRUE,repel=TRUE ,label.size = 0.9))
dev.off()


###harmony predictions
predictions_harmony <- SingleR(test=harmony.sce, assay.type.test=1, 
                              ref=single.sce, labels=single.sce$cell)

table(predictions_harmony$labels)
a <- as.factor(predictions_harmony@listData[["labels"]])
integrated_harmony@meta.data[["pred"]] <- a
b <- integrated_harmony
b@active.ident <- b@meta.data[["pred"]]
saveRDS(b, "./objects/sp/integrated/integrated.harmony_predictions.rds")

pdf(file.path("./results/singleR/",filename = "cell_types_singlecell_harmony.pdf"))
print(SpatialDimPlot(b, combine = FALSE,label=TRUE,repel=TRUE ,label.size = 0.9))
dev.off()


#####Lets check the markers
#Filter
markers_x <- subset(NicheMarkers10x , 0.5 < avg_logFC)

#Top5
top10_x <- markers_x %>%
  group_by(cluster) %>%
  top_n(n = 10,
        wt = avg_logFC)

####Markers per cluster

top5_x$cluster <- as.character(top5_x$cluster)
Schwann_cells <- top5_x[top5_x$cluster %in% c("Schwann cells"), ]
Schwann_cells <- as.vector(Schwann_cells$gene)

cell_types <- as.list(as.vector(unique(NicheMarkers10x$cluster)))
names(cell_types) <- as.vector(unique(NicheMarkers10x$cluster))

for (i in 1:length(cell_types)){
  a <- cell_types[[i]]
  top10_x$cluster <- as.character(top10_x$cluster)
  b <- top10_x[top10_x$cluster %in% c(a), ]
  b <- as.vector(b$gene)
  cell_types[[i]] <- b
}

names(cell_types)[22] <- "Ng2 MSCs"

####AddmoduleScore Seurat

integrated_seurat <- readRDS("./objects/sp/integrated/integrated.seurat_type.rds")
integrated_harmony <- readRDS("./objects/sp/integrated/integrated.harmony_type.rds")

for (i in 1:length(cell_types)){
  a <- cell_types[[i]]
  c <- paste0(names(cell_types[i]), "1")
  c <- gsub("-", ".", c)
  c <- gsub(" ", ".", c)
  c <- gsub("/", ".", c)
  integrated_seurat <- AddModuleScore(integrated_seurat, features = cell_types[i], name = names(cell_types[i]))
  #spatial 
  pdf(file.path("./results/singleR/module_score/",filename = paste0("spatial_seurat_",c,".pdf")))
  print(SpatialFeaturePlot(integrated_seurat,features=c,combine = FALSE))
  dev.off()
}

####AddmoduleScore Harmony

for (i in 1:length(cell_types)){
  a <- cell_types[[i]]
  c <- paste0(names(cell_types[i]), "1")
  c <- gsub("-", ".", c)
  c <- gsub(" ", ".", c)
  c <- gsub("/", ".", c)
  integrated_harmony <- AddModuleScore(integrated_harmony, features = cell_types[i], name = names(cell_types[i]))
  #spatial 
  pdf(file.path("./results/singleR/module_score/",filename = paste0("spatial_harmony_",c,".pdf")))
  print(SpatialFeaturePlot(integrated_harmony,features=c,combine = FALSE))
  dev.off()
}

####Are the different cell types different expressed between conditions?


####AddmoduleScore
integrated_harmony <- AddModuleScore(integrated_harmony, features = list(cell_types[2]), name = names(cell_types[1]))
integrated_seurat <- AddModuleScore(integrated_seurat, features = list(cell_types[2]), name = names(cell_types[1]))

####Plots

pdf(file.path("./results/singleR/module_score/",filename = "Adipo-CAR_singlecell_harmony.pdf"))
print(SpatialFeaturePlot(integrated_harmony,features=c("Adipo.CAR1"),combine = FALSE))
dev.off()

pdf(file.path("./results/singleR/module_score/",filename = "Adipo-CAR_singlecell_seurat.pdf"))
print(SpatialFeaturePlot(integrated_seurat,features=c("Adipo.CAR1"), combine = FALSE))
dev.off()



