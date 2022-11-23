## SCRIPT: Singel cell reference data preparing BM project

## 23.11.22 Laura Sudupe , git @lsudupe

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

a <- as.factor(single_cell_bonemarrow@meta.data[["cell"]])
levels(a) <- list(PSC  = "Ery/Mk prog.", 
                  PSC = "Neutro prog.",
                  PSC = "Mono prog.",
                  PSC = "Gran/Mono prog.",
                  PSC = "LMPPs",
                  Bcells = "large pre-B.",
                  PSC = "Mk prog.",
                  Erythroblasts = "Erythroblasts",
                  PSC = "Eo/Baso prog.",
                  Monocytes = "Monocytes",
                  PSC = "Ery prog.",
                  Bcells = "pro-B",
                  Tcells = "T cells",
                  Neutrophils = "Neutrophils",
                  MSC = "Adipo-CAR",
                  MSC = "Ng2+ MSCs",
                  NC = "Schwann cells",
                  MSC = "Osteoblasts",
                  MSC_fibro = "Arteriolar fibro.",
                  EC = "Sinusoidal ECs",
                  MSC = "Osteo-CAR",
                  Bcell = "small pre-B.",
                  Chondrocytes = "Chondrocytes",
                  MSC_fibro = "Endosteal fibro.",
                  Fibro_Chondro_p = "Fibro/Chondro p.",
                  MSC_fibro = "Stromal fibro.",
                  EC = "Arteriolar ECs",
                  MSC_fibro = "Myofibroblasts",
                  MSC = "Smooth muscle",
                  DC = "Dendritic cells",
                  IC = "NK cells",
                  Bcells = "B cell"
)

#SMC, NC, Fibro_Chondro_p take out, 
single_cell_bonemarrow@meta.data[["ident"]] <- a
Seurat::Idents(object = single_cell_bonemarrow) <- single_cell_bonemarrow@meta.data[["ident"]]

single_cell_bonemarrow <- subset(x = single_cell_bonemarrow, idents = c("NC", "Fibro_Chondro_p"), invert = TRUE)
single_cell_bonemarrow@meta.data[["ident"]] <- single_cell_bonemarrow@active.ident


saveRDS(single_cell_bonemarrow, "./objects/heterogeneity/single_cell_bonemarrow.rds")
single_cell_bonemarrow <- readRDS("./objects/heterogeneity/single_cell_bonemarrow.rds")







