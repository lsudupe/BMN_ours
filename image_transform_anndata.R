## SCRIPT: prepare metadata for python BM project

## 29.03.23 Laura Sudupe , git @lsudupe

## Libraries
library(Seurat)
library(tidyverse)
library(STutility)

##read data
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")
se_healthy <- readRDS("./objects/heterogeneity/healthy/se_hierarchical.rds")

## set ident
se <- SetIdent(se, value = se@meta.data[["name"]])
se_healthy <- SetIdent(se_healthy, value = se_healthy@meta.data[["name"]])
##subset data
M1_fem_1C <- SubsetSTData(se, ident = ("M1_fem_1C"))
M2_F_2B <- SubsetSTData(se, ident = ("M2_F_2B"))
M8_F2_1C <- SubsetSTData(se, ident = ("M8_F2_1C"))
M9_F2_1C <- SubsetSTData(se, ident = ("M9_F2_1C"))
M3_F_1C <- SubsetSTData(se_healthy, ident = ("M3_F_1C"))
M3_fem_1C <- SubsetSTData(se_healthy, ident = ("M3_fem_1C"))


##save data
saveRDS(M1_fem_1C, "./objects/sc/integrated/se_deco_M1_fem_1C.rds")
saveRDS(M2_F_2B,"./objects/sc/integrated/se_deco_M2_F_2B.rds")
saveRDS(M8_F2_1C,"./objects/sc/integrated/se_deco_M8_F2_1C.rds")
saveRDS(M9_F2_1C,"./objects/sc/integrated/se_deco_M9_F2_1C.rds")
saveRDS(M3_F_1C,"./objects/sc/integrated/se_deco_M3_F_1C.rds")
saveRDS(M3_fem_1C,"./objects/sc/integrated/se_deco_M3_fem_1C.rds")

##read data

## Create each individual Seurat object for every sample
samples <- c("M1_fem_1C", "M2_F_2B","M8_F2_1C", "M9_F2_1C", "M3_F_1C" ,  "M3_fem_1C")
lista <- c(M1_fem_1C, M2_F_2B,M8_F2_1C, M9_F2_1C, M3_F_1C , M3_fem_1C)
names(lista) <- samples
prueba <-c()

for (i in 1:length(lista)){
  a <- lista[[i]]
  b <- names(lista[i])
  a_normal <- Load10X_Spatial(data.dir = paste0("./data/data/", names(lista[i])),
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial",
                       slice = names(lista[i]),
                       filter.matrix = TRUE)
  
  ## coordinates
  coor <- a_normal@images[[i]]@coordinates
  ## image data
  staffli_image <- a@tools[["Staffli"]]@meta.data
  staffli_image$barcode <- rownames(staffli_image)
  sample_coor <- staffli_image[3:4]
  
  ##tissue position list
  tissue_position_list <- read.csv(paste0("./data/data/", names(lista[i]),"/spatial/tissue_positions_list.csv"))
  prueba <- tissue_position_list
  coor$barcode <- rownames(coor)
  #coor$barcode <- paste0(coor$barcode, "_5")
  rownames(prueba) <- prueba$barcode
  
  ## order barcodes
  prueba_ <- coor[coor$barcode %in% staffli_image$barcode,]
  
  prueba_$origbarcode <- rownames(prueba_)
  rownames(prueba_) <- prueba_$barcode
  aaa <- c("ee")
  target <- c(rownames(staffli_image))
  prueba_ <- prueba_[match(target, prueba_$barcode),]
  
  ## order data

  rownames(prueba_) <- prueba_$origbarcode
  prueba_$barcode <- NULL
  prueba_$origbarcode <- NULL
  
  ## save data
  write.csv(prueba_, paste0("./data/data/image_info/", names(lista[i]) , "_coordinates.csv"))
  
}



###normal spatial
M3_F_1C_normal <- Load10X_Spatial("./data/data/M3_F_1C/",
                     filename = "filtered_feature_bc_matrix.h5",
                     assay = "Spatial",
                     slice = "M3_F_1C",
                     filter.matrix = TRUE)
coor <- M3_F_1C_normal@images[["M3_F_1C"]]@coordinates

##image data
M3_F_1C_image <- M3_F_1C@tools[["Staffli"]]@meta.data
M3_F_1C_image$barcode <- rownames(M3_F_1C_image)
M3 <- M3_F_1C_image[3:4]
#M8_count <- as.data.frame(M8_F2_1C@assays[["SCT"]]@data)

##tissue position list
tissue_position_list_M3 <- read.csv("./data/data/M3_F_1C/spatial/tissue_positions_list.csv")
prueba <- tissue_position_list_M3
coor$barcode <- rownames(coor)
coor$barcode <- paste0(coor$barcode, "_3")
rownames(prueba) <- prueba$barcode
#prueba$barcode <- NULL
prueba_ <- coor[coor$barcode %in% M3_F_1C_image$barcode,]
#prueba_ <- prueba_ %>% mutate(barcode = gsub("_5", "", barcode))
##order barcodes
prueba_$origbarcode <- rownames(prueba_)
rownames(prueba_) <- prueba_$barcode
target <- c(rownames(M3_F_1C_image))
prueba_ <- prueba_[match(target, prueba_$barcode),]
##
rownames(prueba_) <- prueba_$origbarcode
prueba_$barcode <- NULL
prueba_$origbarcode <- NULL

write.csv(prueba_, "./data/data/image_info/M3_coordinates.csv")



write.csv(M8_count, "data/data/image_info/M8_count.csv")
write.csv(M8, "data/data/image_info/M8.csv")








