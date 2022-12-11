## SCRIPT: Clustering deconvolution results BM project

## 11.12.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)

#Data---------------------------------
femur <- readRDS("./objects/card/heterogeneity/M1_fem_1C_subgroup_deco.rds")


#Analysis---------------------------------
meta <- femur@meta.data
types <- meta[,11:22]

