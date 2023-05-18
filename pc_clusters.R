## SCRIPT: PC cluster differentiation visium data BM project

## 18.05.23 Laura Sudupe , git @lsudupe

## Libraries
library(Seurat)
library(tidyverse)
library(STutility)

# Data
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

Idents(object = se) <- "name"
#M8
M8 <- SubsetSTData(se, ident = "M8_F2_1C")
meta <- M8@meta.data
rownames(meta) <- sub("_5", "", rownames(meta))
M8@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m8_clustering.csv", row.names=TRUE)
M8_s <- ManualAnnotation(M8)
#M2
M2 <- SubsetSTData(se, ident = "M2_F_2B")
meta <- M2@meta.data
rownames(meta) <- sub("_5", "", rownames(meta))
M2@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m2_clustering.csv", row.names=TRUE)
M2_s <- ManualAnnotation(M2)
