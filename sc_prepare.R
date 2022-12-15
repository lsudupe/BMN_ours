## SCRIPT: Singel cell reference data preparing BM project

## 23.11.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(scCustomize)
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

a <- as.factor(single_cell_bonemarrow@meta.data[["cell"]])
levels(a) <- list(PSC  = "Ery/Mk prog.", 
                  PSC = "Neutro prog.",
                  PSC = "Mono prog.",
                  PSC = "Gran/Mono prog.",
                  PSC = "LMPPs",
                  Bcell = "large pre-B.",
                  PSC = "Mk prog.",
                  Erythroblasts = "Erythroblasts",
                  PSC = "Eo/Baso prog.",
                  Monocytes = "Monocytes",
                  PSC = "Ery prog.",
                  Bcell = "pro-B",
                  Tcell = "T cells",
                  Neutrophils = "Neutrophils",
                  MSC = "Adipo-CAR",
                  MSC = "Ng2+ MSCs",
                  NC = "Schwann cells",
                  MSC = "Osteoblasts",
                  Fibro = "Arteriolar fibro.",
                  EC = "Sinusoidal ECs",
                  MSC = "Osteo-CAR",
                  Bcell = "small pre-B.",
                  Chondrocytes = "Chondrocytes",
                  Fibro = "Endosteal fibro.",
                  Fibro = "Fibro/Chondro p.",
                  Fibro = "Stromal fibro.",
                  EC = "Arteriolar ECs",
                  Fibro = "Myofibroblasts",
                  MSC = "Smooth muscle",
                  DC = "Dendritic cells",
                  NK = "NK cells",
                  Bcell = "B cell"
)

#SMC, NC, Fibro_Chondro_p take out, 
single_cell_bonemarrow@meta.data[["ident"]] <- a
Seurat::Idents(object = single_cell_bonemarrow) <- single_cell_bonemarrow@meta.data[["ident"]]

single_cell_bonemarrow <- subset(x = single_cell_bonemarrow, idents = c("NC", "Fibro_Chondro_p"), invert = TRUE)
single_cell_bonemarrow@meta.data[["ident"]] <- single_cell_bonemarrow@active.ident


saveRDS(single_cell_bonemarrow, "./objects/heterogeneity/single_cell_bonemarrow.rds")
single_cell_bonemarrow <- readRDS("./objects/heterogeneity/single_cell_bonemarrow.rds")


############READ AZARI data
PC_MM <- readRDS("./data/single-cell/PC/scRNA_MM_PC.rds")
PC_MM <- subset(x = PC_MM, idents = c("MM_MIC"))
PC_MM@meta.data[["ident"]] <- PC_MM@active.ident


#####integrate data
##Merge them
combined_sc <- merge(single_cell_bonemarrow, y = c( PC_MM), 
                  add.cell.ids = c("single_cell_bonemarrow", "PC_MM"), project = "BM")
combined_sc@meta.data[["split"]] <- combined_sc@active.ident
x <- combined_sc
# Select the most variable features to use for integration

list <- SplitObject(x, split.by = "split")
list <- lapply(X = list, FUN = SCTransform, assay="RNA")
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
x <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

x <- RunPCA(x, verbose=FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20) #%>%
  #FindNeighbors(reduction = "pca", dims = 1:20) %>%
  #FindClusters(resolution=0.8) 

####marker genes for MM
# Determine differentiating markers for PC_MM
Seurat::Idents(object = PC_MM) <- PC_MM@meta.data[["orig.ident"]]
PC_MM <- PrepSCTFindMarkers(PC_MM, assay = "SCT", verbose = TRUE)
MM_genes <- Seurat::FindAllMarkers(object = PC_MM, 
                                              assay = "SCT",
                                              verbose = TRUE, 
                                              only.pos = TRUE)
cluster_MIC <- subset(MM_genes, p_val_adj < 0.05 & 0.05 < avg_log2FC)
genes <- cluster_MIC[grepl("MM_MIC", cluster_MIC[,6]),]
genes <- cluster_MIC$gene

#Top5
top50 <- cluster_6_f %>%
  top_n(n = 50,
        wt = avg_log2FC)

write.csv(genes, "./data/single-cell/PC/MM_MIC_genes.csv", row.names = FALSE)




