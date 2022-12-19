## SCRIPT: Jackard index deconvolution hyerarquical and seurat clusters femur BM project

## 19.12.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2) 
library(scclusteval)


#Data--------------------------------------
integrated_seurat <- readRDS("./objects/heterogeneity/femur_hierarchical.rds") #femur
x <- integrated_seurat
x <- SetIdent(x, value = x@meta.data[["seurat_clusters"]])
z <- integrated_seurat
z <- SetIdent(z, value = z@meta.data[["clustering"]])




######Jackard index 
pdf(file.path("./results/jackard/femur/",filename = "jackard_femur.pdf"))                                                         
PairWiseJaccardSetsHeatmap(x@active.ident,
                           z@active.ident,
                           best_match = TRUE,
                           col_low = "white",
                           col_high = "red"
                           )
dev.off()
