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
#M8_s <- ManualAnnotation(M8)
meta <- M8_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M8_", labels, "_cluster", clustering))
M8_s@meta.data <- meta
saveRDS(M8_s, "./objects/pc_clusters/M8_s.rds")

#M2
M2 <- SubsetSTData(se, ident = "M2_F_2B")
meta <- M2@meta.data
rownames(meta) <- sub("_5", "", rownames(meta))
M2@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m2_clustering.csv", row.names=TRUE)
#M2_s <- ManualAnnotation(M2)
meta <- M2_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M2_", labels, "_cluster", clustering))
M2_s@meta.data <- meta
saveRDS(M2_s, "./objects/pc_clusters/M2_s.rds")

#M9
M9 <- SubsetSTData(se, ident = "M9_F2_1C")
meta <- M9@meta.data
rownames(meta) <- sub("_5", "", rownames(meta))
M9@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m9_clustering.csv", row.names=TRUE)
M9_s <- ManualAnnotation(M9)
meta <- M9_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M9_", labels, "_cluster", clustering))
M9_s@meta.data <- meta
saveRDS(M9_s, "./objects/pc_clusters/M9_s.rds")

#M1
M1 <- SubsetSTData(se, ident = "M1_fem_1C")
meta <- M1@meta.data
rownames(meta) <- sub("_5", "", rownames(meta))
M1@meta.data <- meta
meta <- meta[, ncol(meta), drop = FALSE]

write.csv(meta, "./data/areas/m1_clustering.csv", row.names=TRUE)
#M1_s <- ManualAnnotation(M1)
meta <- M1_s@meta.data
meta <- meta %>%
  mutate(pc_clusters = paste0("M1_", labels, "_cluster", clustering))
M1_s@meta.data <- meta
saveRDS(M1_s, "./objects/pc_clusters/M1_s.rds")

###combine the data

se_s <- MergeSTData(M1_s, y = c(M2_s, M8_s, M9_s), 
                  add.spot.ids = c("M1_s", "M2_s", "M8_s", "M9_s"), project = "BM")

saveRDS(se_s, "./objects/pc_clusters/combined_s.rds")

################ CHECK DORMANT SIGNATURE





