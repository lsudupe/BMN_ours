## SCRIPT: Cell type identification Itziar canonical markers in the spatial data BMN project

## 25.09.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2) 
library('magrittr')
library(tidyverse)
library(base)

###Data
canonical <- read.csv("./data/markers_canonical.csv",sep=";")

integrated_seurat <- readRDS("./objects/sp/integrated/integrated.seurat_type.rds")
integrated_harmony <- readRDS("./objects/sp/integrated/integrated.harmony_type.rds")

####AddmoduleScore Seurat
lis <- canonical
lis$RIBOSOMAL <- NULL
lis <- as.list(lis)

for (i in 1:length(lis)){
  a <- lis[[i]]
  c <- paste0(names(lis[i]), "1")
  integrated_seurat <- AddModuleScore(integrated_seurat, features = lis[i], name = names(lis[i]))
  #spatial 
  pdf(file.path("./results/marker_genes/Itziar_canonical/seurat/",filename = paste0("spatial_seurat_",c,".pdf")))
  print(SpatialFeaturePlot(integrated_seurat,features=c,combine = FALSE))
  dev.off()
}

####AddmoduleScore Harmony
canonical <- as.list(canonical)
for (i in 1:length(canonical)){
  a <- cell_types[[i]]
  c <- paste0(names(canonical[i]), "1")
  integrated_harmony <- AddModuleScore(integrated_harmony, features = canonical[i], name = names(canonical[i]))
  #spatial 
  pdf(file.path("./results/marker_genes/Itziar_canonical/harmony/",filename = paste0("spatial_harmony_",c,".pdf")))
  print(SpatialFeaturePlot(integrated_harmony,features=c,combine = FALSE))
  dev.off()
}

