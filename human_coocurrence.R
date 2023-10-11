## SCRIPT: Human sample pc, tex, reads co-ocurrence

## 11.10.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(Signac)
library(RColorBrewer)
library(STutility)
library(UCell)
library(readr)
library(dplyr)
library(gridExtra)
library(tibble)

#pc_markers
pc_mm_markers <- c("SHH", "DHH", "IHH", "PTCH1", "PTCH2", "SMO", "SUFU", "GLI1", 
                   "GLI2", "GLI3", "CD19", "CD44", "CXCR4", "KLF4", "CD28", "CD33", "CD27")

human_tcell_exhausted <- c("PDCD1", "CTLA4", "TNFRSF9", "HAVCR2", "TOX", "TIGIT", "WARS", "RSAD2",
                           "MCM7", "MX1", "NDFIP2", "ENOSF1", "CCDC141", "STMN1", "TTN", "FASLG",
                           "MCM5", "NAB1", "PHLDA1", "MCM3", "PCNA", "GAPDH", "OASL", "IFI44L",
                           "TBC1D4", "SLC43A3", "PAM", "CCL3", "ACP5", "OAS3", "CD38", "TNFSF10",
                           "GBP2", "KIF20B", "CTSB")

human <- readRDS("./objects/sp/human/human_combined.rds")


# split human
B08041 <- human[["BM_human_AP-B08041_"]]
B08805 <- human[["BM_human_AP-B08805"]]
B10395 <- human[["BM_B10395"]]

lista <- c(B08041, B08805, B10395)
names(lista) <- c("B08041","B08805", "B10395")

###Enrichment score
color <- brewer.pal(11,"Spectral")
color <- rev(color)

post <- c()
for (i in 1:length(lista)){
  a <- lista[[i]]

  # Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(pc_mm_markers))
  a@meta.data[["signature_1_pc"]] <- as.vector(vector)
    
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(human_tcell_exhausted))
  a@meta.data[["signature_1_tex"]] <- as.vector(vector)
  
  post[[length(post) + 1]] <- a
}

names(post) <- c("B08041","B08805", "B10395")


meta <- post[["B08041"]]@meta.data
