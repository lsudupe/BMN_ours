## SCRIPT: Spatial object cluster comparation BM project

## 22.09.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(scCustomize)
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

pdf(file.path("./results/clusters/single_cell/",filename = "cell_types_singlecell.pdf"))
#print(DimPlot(single_cell_bonemarrow, reduction = "tsne", label = TRUE))
print(DimPlot_scCustom(single_cell_bonemarrow, reduction = "tsne", label = TRUE, repel=TRUE))
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
saveRDS(predictions_seurat, "./objects/sp/integrated/seurat_predictions.rds")

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
saveRDS(predictions_harmony, "./objects/sp/integrated/harmony_predictions.rds")

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

seurat <- readRDS("./objects/sp/integrated/integrated.seurat_type.rds")
harmony <- readRDS("./objects/sp/integrated/integrated.harmony_type.rds")

seurat_cell_types <- cell_types
seurat_cell_types[["T cells"]] <- NULL

for (i in 1:length(seurat_cell_types)){
  a <- seurat_cell_types[[i]]
  c <- paste0(names(seurat_cell_types[i]), "1")
  c <- gsub("-", ".", c)
  c <- gsub(" ", ".", c)
  c <- gsub("/", ".", c)
  seurat <- AddModuleScore(seurat, features = seurat_cell_types[i], name = names(seurat_cell_types[i]))
  #spatial 
  pdf(file.path("./results/singleR/module_score/",filename = paste0("spatial_seurat_",c,".pdf")))
  print(SpatialFeaturePlot(seurat,features=c,combine = FALSE))
  dev.off()
}

####AddmoduleScore Harmony
harmony_cell_types <- cell_types

for (i in 1:length(harmony_cell_types)){
  a <- harmony_cell_types[[i]]
  c <- paste0(names(harmony_cell_types[i]), "1")
  c <- gsub("-", ".", c)
  c <- gsub(" ", ".", c)
  c <- gsub("/", ".", c)
  harmony <- AddModuleScore(harmony, features = harmony_cell_types[i], name = names(harmony_cell_types[i]))
  #spatial 
  pdf(file.path("./results/singleR/module_score/",filename = paste0("spatial_harmony_",c,".pdf")))
  print(SpatialFeaturePlot(harmony,features=c,combine = FALSE))
  dev.off()
}

####Are the different cell types different expressed between conditions?

pdf(file.path("./results/singleR/",filename = "harmony_heatmap.pdf"))
plotScoreHeatmap(predictions_harmony)
dev.off()

pdf(file.path("./results/singleR/",filename = "harmony_violin.pdf"))
plotDeltaDistribution(predictions_harmony, ncol = 3)
dev.off()

pdf(file.path("./results/singleR/",filename = "seurat_heatmap.pdf"))
plotScoreHeatmap(predictions_seurat)
dev.off()

pdf(file.path("./results/singleR/",filename = "seurat_violin.pdf"))
plotDeltaDistribution(predictions_seurat, ncol = 3)
dev.off()

###markers
library(scater)
all.markers_seurat <- metadata(predictions_seurat)$de.genes
seurat.sce$labels <- predictions_seurat$labels

pdf(file.path("./results/singleR/",filename = "seurat_heatmap.pdf"))
plotHeatmap(seurat.sce, order_columns_by="labels",
            features=unique(unlist(all.markers_seurat$Chondrocytes)),
            colour_columns_by=c("MM", "control")) 
dev.off()


all.markers_harmony <- metadata(predictions_harmony)$de.genes
harmony.sce$labels <- predictions_harmony$labels

pdf(file.path("./results/singleR/",filename = "harmony_heatmap.pdf"))
plotHeatmap(harmony.sce, order_columns_by="labels",
            features=unique(unlist(all.markers_harmony$Chondrocytes)))
dev.off()

####SCORES
scores_seurat <- as.data.frame(predictions_seurat@listData[["scores"]])

first_max <- c()
second_max <- c()
rest_seurat <- c()
for (i in 1:nrow(scores_seurat)) {
  first_max = c(first_max, max(scores_seurat[i,]))
  second_max = c(second_max, max(scores_seurat[i,][scores_seurat[i,] != max(scores_seurat[i,])]))
  rest_seurat <- (first_max - second_max)
}

library(hrbrthemes)
integrated_seurat@meta.data[["SingleR_score"]] <- first_max

a <-integrated_seurat@meta.data %>% 
  ggplot(aes(x= SingleR_score, y=nFeature_Spatial,color=pred,  label=pred)) + 
  geom_point(size=0.6) 
b <-integrated_seurat@meta.data %>% 
  ggplot(aes(x= SingleR_score, y=nCount_Spatial,color=pred, label=pred)) + 
  geom_point(size=0.6) 


pdf(file.path("./results/singleR/",filename = "seurat_score_feature.pdf"))
plot(a)
dev.off()

pdf(file.path("./results/singleR/",filename = "seurat_score_umi.pdf"))
plot(b)
dev.off()

