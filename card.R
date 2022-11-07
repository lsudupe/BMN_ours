## SCRIPT: Deconvolution using CARD areas BM project

## 30.10.22 Laura Sudupe , git @lsudupe
#https://github.com/YingMa0107/CARD/

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
#devtools::install_github('YingMa0107/CARD')
library(CARD)
library(gtools)
library(scatterpie)
library(ggcorrplot)

source(file = "./card.plot2.R")


#Data---------------------------------
spatial <- readRDS("./objects/sp/second/combined_filtered.rds")
single_cell_bonemarrow <- readRDS("./objects/sc/single_cell_bonemarrow.rds")


######Prepare data
#spatial
x <- spatial
x.image <- x@images
x@images[["M1_tib_1A"]]<- NULL
x@images[["M1_fem_1C"]]<- NULL
x@images[["M3_tib_2A"]]<- NULL
x@images[["M3_fem_1C"]]<- NULL
list <- SplitObject(x, split.by = "orig.ident")

#sc
single_cell_bonemarrow_counts <- single_cell_bonemarrow@assays[["RNA"]]@counts
single_cell_bonemarrow_meta <- single_cell_bonemarrow@meta.data
#single_cell_bonemarrow_meta <- single_cell_bonemarrow_meta[, 6:7]

for (i in 1:length(list)){
  a <- list[[i]]
  a <- SCTransform(a, assay="Spatial")
  a <- RunPCA(a, npcs = 20, verbose = FALSE, assay="SCT")
  a <- RunUMAP(a, reduction = "pca", dims = 1:20, verbose = FALSE)
  a <- FindNeighbors(a, reduction = "pca", dims = 1:20)
  a <- FindClusters(a, resolution=0.5)
  a_count <- a@assays[["Spatial"]]@counts
  a_location <- x.image[[names(list[i])]]@coordinates
  a_location <- a_location[,2:3]
  colnames(a_location) <- c("x", "y")
  ###################CARD object creation
  CARD_obj_a = createCARDObject(
    sc_count = single_cell_bonemarrow_counts,
    sc_meta = single_cell_bonemarrow_meta,
    spatial_count = a_count,
    spatial_location = a_location,
    ct.varname = "groups",
    ct.select = unique(single_cell_bonemarrow_meta$groups),
    sample.varname = "metadata....experiment..",
    minCountGene = 100,
    minCountSpot = 5)
  ###################Deco
  CARD_a = CARD_deconvolution(CARD_object = CARD_obj_a)
  CARD_a_proportions <- CARD_a@Proportion_CARD
  ###################Plots
  p1 <- CARD.visualize.pie(proportion = CARD_a@Proportion_CARD,
                           spatial_location = CARD_a@spatial_location)
  pdf(paste("./results/CARD/", names(list[i]),"_1.pdf",sep=""))
  print(p1)
  dev.off()
  ## select the cell type that we are interested
  ct.visualize = c("HSC_PSC","IC","MSC","NC","EC")
  ## visualize the spatial distribution of the cell type proportion
  p2 <- CARD.visualize.prop.2(
    proportion = CARD_a@Proportion_CARD,        
    spatial_location = CARD_a@spatial_location, 
    ct.visualize = ct.visualize,                 ### selected cell types to visualize
    colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
    NumCols = 4)                                 ### number of columns in the figure panel
  pdf(paste("./results/CARD/", names(list[i]),"_2.pdf",sep=""))
  print(p2)
  dev.off()
  ## correlation
  p3 <- CARD.visualize.Cor(CARD_a@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
  pdf(paste("./results/CARD/", names(list[i]),"_3.pdf",sep=""))
  print(p3)
  dev.off()
  ## save object
  saveRDS(a,file = paste0("./objects/card/",names(list[i]),".rds"))
  ## save card results
  saveRDS(CARD_a,file = paste0("./objects/card/",names(list[i]),"_CARD_obj.rds"))
}

## Card object
M1_fem_1C_CARD_obj <- readRDS("./objects/card/M1_fem_1C_CARD_obj.rds")
M1_tib_1A_CARD_obj <- readRDS("./objects/card/M1_tib_1A_CARD_obj.rds")
M3_fem_1C_CARD_obj <- readRDS("./objects/card/M3_fem_1C_CARD_obj.rds")
M3_tib_2A_CARD_obj <- readRDS("./objects/card/M3_tib_2A_CARD_obj.rds")

card_list <- c(M1_fem_1C_CARD_obj, M1_tib_1A_CARD_obj, M3_fem_1C_CARD_obj, M3_tib_2A_CARD_obj)
names(card_list) <- c("M1_fem_1C_CARD_obj","M1_tib_1A_CARD_obj","M3_fem_1C_CARD_obj","M3_tib_2A_CARD_obj")

card_list <- c(M1_fem_1C_CARD_obj, M1_tib_1A_CARD_obj)
names(card_list) <- c("M1_fem_1C_CARD_obj","M1_tib_1A_CARD_obj")

## Spatial object
M1_fem_1C <- readRDS("./objects/card/M1_fem_1C.rds")
M1_tib_1A <- readRDS("./objects/card/M1_tib_1A.rds")
M3_fem_1C <- readRDS("./objects/card/M3_fem_1C.rds")
M3_tib_2A <- readRDS("./objects/card/M3_tib_2A.rds")

spatial_list <- c(M1_fem_1C, M1_tib_1A, M3_fem_1C, M3_tib_2A)
names(spatial_list) <- c("M1_fem_1C","M1_tib_1A","M3_fem_1C","M3_tib_2A")

spatial_list <- c(M1_fem_1C, M1_tib_1A)
names(spatial_list) <- c("M1_fem_1C","M1_tib_1A")

list <- list(spatial_list, card_list)
names(list) <- c("spatial_list", "card_list")
  
######Proportions#######################################################################
#https://stackoverflow.com/questions/46205479/looping-over-multiple-lists-with-base-r
for (i in 1:length(list)){
  a <- list[[1]][[i]]
  cc <- list[[1]][i]
  c <- names(cc)
  b <- list[[2]][[i]]
  ## 
  a_meta <- a@meta.data
  a_area <- a_meta$area
  a_proportions <- b@Proportion_CARD
  a_proportions <- as.data.frame(a_proportions)
  a_proportions["area"] <- as.vector(a_area)
  ## proportions per cell type and area
  # AC
  AC_value <- a_proportions[grepl("AC", a_proportions[,6]),]
  AC_value$area <- NULL
  AC_HSC_PSC <- (sum(AC_value$HSC_PSC))*100/nrow(AC_value)
  AC_MSC <- (sum(AC_value$MSC))*100/nrow(AC_value)
  AC_NC <- (sum(AC_value$NC))*100/nrow(AC_value)
  AC_EC <- (sum(AC_value$EC))*100/nrow(AC_value)
  AC_IC <- (sum(AC_value$IC))*100/nrow(AC_value)
  AC_proportions <- c(AC_HSC_PSC, AC_MSC , AC_NC, AC_EC, AC_IC)
  # BM
  BM_value <- a_proportions[grepl("BM", a_proportions[,6]),]
  BM_value$area <- NULL
  BM_HSC_PSC <- (sum(BM_value$HSC_PSC))*100/nrow(BM_value)
  BM_MSC <- (sum(BM_value$MSC))*100/nrow(BM_value)
  BM_NC <- (sum(BM_value$NC))*100/nrow(BM_value)
  BM_EC <- (sum(BM_value$EC))*100/nrow(BM_value)
  BM_IC <- (sum(BM_value$IC))*100/nrow(BM_value)
  BM_proportions <- c(BM_HSC_PSC, BM_MSC , BM_NC, BM_EC, BM_IC)
  # Muscle
  Muscle_value <- a_proportions[grepl("Muscle", a_proportions[,6]),]
  Muscle_value$area <- NULL
  Muscle_HSC_PSC <- (sum(Muscle_value$HSC_PSC))*100/nrow(Muscle_value)
  Muscle_MSC <- (sum(Muscle_value$MSC))*100/nrow(Muscle_value)
  Muscle_NC <- (sum(Muscle_value$NC))*100/nrow(Muscle_value)
  Muscle_EC <- (sum(Muscle_value$EC))*100/nrow(Muscle_value)
  Muscle_IC <- (sum(Muscle_value$IC))*100/nrow(Muscle_value)
  Muscle_proportions <- c(Muscle_HSC_PSC, Muscle_MSC , Muscle_NC, Muscle_EC, Muscle_IC)
  #GP
  GP_value <- a_proportions[grepl("GP", a_proportions[,6]),]
  GP_value$area <- NULL
  GP_HSC_PSC <- (sum(GP_value$HSC_PSC))*100/nrow(GP_value)
  GP_MSC <- (sum(GP_value$MSC))*100/nrow(GP_value)
  GP_NC <- (sum(GP_value$NC))*100/nrow(GP_value)
  GP_EC <- (sum(GP_value$EC))*100/nrow(GP_value)
  GP_IC <- (sum(GP_value$IC))*100/nrow(GP_value)
  GP_proportions <- c(GP_HSC_PSC, GP_MSC , GP_NC, GP_EC, GP_IC)
  # Bone
  Bone_value <- a_proportions[grepl("Bone", a_proportions[,6]),]
  Bone_value$area <- NULL
  Bone_HSC_PSC <- (sum(Bone_value$HSC_PSC))*100/nrow(Bone_value)
  Bone_MSC <- (sum(Bone_value$MSC))*100/nrow(Bone_value)
  Bone_NC <- (sum(Bone_value$NC))*100/nrow(Bone_value)
  Bone_EC <- (sum(Bone_value$EC))*100/nrow(Bone_value)
  Bone_IC <- (sum(Bone_value$IC))*100/nrow(Bone_value)
  Bone_proportions <- c(Bone_HSC_PSC, Bone_MSC , Bone_NC, Bone_EC, Bone_IC)
  ## dataframe
  pro_df_ <- data.frame(AC_proportions, BM_proportions, Muscle_proportions ,GP_proportions, Bone_proportions)
  rownames(pro_df_) <- c("HSC_PSC", "MSC", "NC", "EC", "IC")
  write.csv(pro_df_, file = paste0("./results/CARD/",c,".csv"))
}


M1_tib_1A_meta <- M1_tib_1A@meta.data
M1_tib_1A_area <- M1_tib_1A_meta$area
CARD_M1_tib_1A_proportions <- as.data.frame(CARD_M1_tib_1A_proportions)
CARD_M1_tib_1A_proportions["area"] <- as.vector(M1_tib_1A_area)

########create a function to get proportions

AC_value <- CARD_M1_tib_1A_proportions[grepl("AC", CARD_M1_tib_1A_proportions[,6]),]
AC_value$area <- NULL
AC_HSC_PSC <- (sum(AC_value$HSC_PSC))*100/nrow(AC_value)
AC_MSC <- (sum(AC_value$MSC))*100/nrow(AC_value)
AC_NC <- (sum(AC_value$NC))*100/nrow(AC_value)
AC_EC <- (sum(AC_value$EC))*100/nrow(AC_value)
AC_EC <- (sum(AC_value$EC))*100/nrow(AC_value)
AC_proportions <- c(AC_HSC_PSC, AC_MSC , AC_NC, AC_EC)

BM_value <- CARD_M1_tib_1A_proportions[grepl("BM", CARD_M1_tib_1A_proportions[,6]),]
BM_value$area <- NULL
BM_HSC_PSC <- (sum(BM_value$HSC_PSC))*100/nrow(BM_value)
BM_MSC <- (sum(BM_value$MSC))*100/nrow(BM_value)
BM_NC <- (sum(BM_value$NC))*100/nrow(BM_value)
BM_EC <- (sum(BM_value$EC))*100/nrow(BM_value)
BM_proportions <- c(BM_HSC_PSC, BM_MSC , BM_NC, BM_EC)

Muscle_value <- CARD_M1_tib_1A_proportions[grepl("Muscle", CARD_M1_tib_1A_proportions[,6]),]
Muscle_value$area <- NULL
Muscle_HSC_PSC <- (sum(Muscle_value$HSC_PSC))*100/nrow(Muscle_value)
Muscle_MSC <- (sum(Muscle_value$MSC))*100/nrow(Muscle_value)
Muscle_NC <- (sum(Muscle_value$NC))*100/nrow(Muscle_value)
Muscle_EC <- (sum(Muscle_value$EC))*100/nrow(Muscle_value)
Muscle_proportions <- c(Muscle_HSC_PSC, Muscle_MSC , Muscle_NC, Muscle_EC)

GP_value <- CARD_M1_tib_1A_proportions[grepl("GP", CARD_M1_tib_1A_proportions[,6]),]
GP_value$area <- NULL
GP_HSC_PSC <- (sum(GP_value$HSC_PSC))*100/nrow(GP_value)
GP_MSC <- (sum(GP_value$MSC))*100/nrow(GP_value)
GP_NC <- (sum(GP_value$NC))*100/nrow(GP_value)
GP_EC <- (sum(GP_value$EC))*100/nrow(GP_value)
GP_proportions <- c(GP_HSC_PSC, GP_MSC , GP_NC, GP_EC)

Bone_value <- CARD_M1_tib_1A_proportions[grepl("Bone", CARD_M1_tib_1A_proportions[,6]),]
Bone_value$area <- NULL
Bone_HSC_PSC <- (sum(Bone_value$HSC_PSC))*100/nrow(Bone_value)
Bone_MSC <- (sum(Bone_value$MSC))*100/nrow(Bone_value)
Bone_NC <- (sum(Bone_value$NC))*100/nrow(Bone_value)
Bone_EC <- (sum(Bone_value$EC))*100/nrow(Bone_value)
Bone_proportions <- c(Bone_HSC_PSC, Bone_MSC , Bone_NC, Bone_EC)

pro_df <- data.frame(AC_proportions, BM_proportions, Muscle_proportions ,GP_proportions, Bone_proportions)
rownames(pro_df) <- c("HSC_PSC", "MSC", "NC", "EC")

as.numeric(rownames(AC_value)) 
nrow(AC_value)

BM_value <- CARD_M1_tib_1A_proportions[grepl("BM", CARD_M1_tib_1A_proportions[,6]),]
BM_value$area <- NULL
Muscle_value <- CARD_M1_tib_1A_proportions[grepl("Muscle", CARD_M1_tib_1A_proportions[,6]),]
Muscle_value$area <- NULL
GP_value <- CARD_M1_tib_1A_proportions[grepl("GP", CARD_M1_tib_1A_proportions[,6]),]
GP_value$area <- NULL
Bone_value <- CARD_M1_tib_1A_proportions[grepl("Bone", CARD_M1_tib_1A_proportions[,6]),]
Bone_value$area <- NULL







