## SCRIPT: Cell type identification in the spatial data BMN project

## 15.09.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2) 
library('magrittr')
library(tidyverse)
library(base)

#Data--------------------------------------
M1_tib_1A  <- readRDS("./objects/sp/M1_tib_1A_normalized.rds")
x <- M1_tib_1A


MSC <- c("Lepr","Cxcl12", "Cd271", "Stro-1","Cd146","Cd90","Thy1", "Lngfr","Cd106","Vcam-1", "Stro-1","Cs146", "Mcam", 
         "Cd105", "Endoglin", "Cd51", "Itgav", "Cd140a", "Pdgfra", "Susd2", "Ly6a", "Cd44", "Ng2", "Nestin", "Kitl", "Nes",
         "Cspg4" )
BMEC <- c("Cd45","Cd235","Cd31","Cd9", "Cdh5","Ly6a", "Pecam1", "Kdr", "Emcn" )
HSC <- c("Cd34", "Sca-1", "Cd27","Cd43","Cd48","Cd117","Cd150", "Cd38", "Cd45RA", "Cd133")
OLCs <- c("Bglap" )
chondrocytes <- c("Acan", "Col2a1", "Sox9", "Runx2", "Ihh", "Mef2c", "Col10a1", "Nt5e", "Cspg4", "Cilp")
fibroblasts <- c("S100a4", "Dcn", "Sema3c", "Cxcl12", "Angpt1", "Sox9")
pericytes <- c("Acta2", "Myh11", "Mcam", "Jag1")

list_top <- list(MSC, BMEC, HSC, OLCs, chondrocytes, fibroblasts, pericytes)
names(list_top) <- c("MSC", "BMEC", "HSC","OLCs", "chondrocytes", "fibroblasts", "pericytes")

for (i in 1:length(list_top)){
  a <- list_top[[i]]
  p1 <- SpatialFeaturePlot(x, features = a, ncol = 3, combine = FALSE)
  pdf(file.path("./results/features/",filename = paste0(names(list_top[i]),".spatial_genes.pdf",sep="")))
  print(p1)
  dev.off()
  
}

EC_arteri <- c("Ly6a", "Ly6c1", "Igfbp3", "Vim", "IgfBp7", "Ppia")
EC_sinusoid <- c("Adamts5", "Stab2", "Il6st", "Ubd", "IgfBp7", "Ppia", "Cd164", "Birv")

list_EC <- list(EC_arteri, EC_sinusoid)
names(list_EC) <- c("EC_arteri", "EC_sinusoid")
             
for (i in 1:length(list_EC)){
  a <- list_EC[[i]]
  p1 <- SpatialFeaturePlot(x, features = a, ncol = 3, combine = FALSE)
  pdf(file.path("./results/features/",filename = paste0(names(list_EC[i]),".spatial_genes.pdf",sep="")))
  print(p1)
  dev.off()
  
}

MSC_jin <- c("Cxcl12", "Lepr", "Igtb1", "Adipoq", "Vcam1")
ONL_jin <- c("Bglap", "Sbds", "Cd200", "Alpl", "Vkorc1", "Col1a1", "Enpp1")

list_MSC <- list(MSC_jin, ONL_jin)
names(list_MSC) <- c("MSC_jin", "ONL_jin")

for (i in 1:length(list_MSC)){
  a <- list_MSC[[i]]
  p1 <- SpatialFeaturePlot(x, features = a, ncol = 3, combine = FALSE)
  pdf(file.path("./results/features/",filename = paste0(names(list_MSC[i]),".spatial_genes.pdf",sep="")))
  print(p1)
  dev.off()
  
}


