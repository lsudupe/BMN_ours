## SCRIPT: prepare metadata for python BM project

## 29.03.23 Laura Sudupe , git @lsudupe

## Libraries
library(Seurat)
library(tidyverse)
library(STutility)

##read data
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")
## set ident
se <- SetIdent(se, value = se@meta.data[["name"]])
##subset data
M1_fem_1C <- SubsetSTData(se, ident = ("M1_fem_1C"))
M2_F_2B <- SubsetSTData(se, ident = ("M2_F_2B"))
M8_F2_1C <- SubsetSTData(se, ident = ("M8_F2_1C"))
M9_F2_1C <- SubsetSTData(se, ident = ("M9_F2_1C"))
##save data
saveRDS(M1_fem_1C, "./objects/sc/integrated/se_deco_M1_fem_1C.rds")
saveRDS(M2_F_2B,"./objects/sc/integrated/se_deco_M2_F_2B.rds")
saveRDS(M8_F2_1C,"./objects/sc/integrated/se_deco_M8_F2_1C.rds")
saveRDS(M9_F2_1C,"./objects/sc/integrated/se_deco_M9_F2_1C.rds")

###normal spatial
M8_F2_1C_normal <- Load10X_Spatial("./data/data/M8_F2_1C/",
                     filename = "filtered_feature_bc_matrix.h5",
                     assay = "Spatial",
                     slice = "M8_F2_1C",
                     filter.matrix = TRUE)
coor <- M8_F2_1C_normal@images[["M8_F2_1C"]]@coordinates


##image data
M8_F2_1C_image <- M8_F2_1C@tools[["Staffli"]]@meta.data
M8_F2_1C_image$barcode <- rownames(M8_F2_1C_image)
M8 <- M8_F2_1C_image[3:4]
#M8_count <- as.data.frame(M8_F2_1C@assays[["SCT"]]@data)

##tissue position list
tissue_position_list_M8 <- read.csv("./data/data/M8_F2_1C/spatial/tissue_positions_list.csv")
prueba <- tissue_position_list_M8
coor$barcode <- rownames(coor)
coor$barcode <- paste0(coor$barcode, "_5")
rownames(prueba) <- prueba$barcode
#prueba$barcode <- NULL
prueba_ <- coor[coor$barcode %in% M8_F2_1C_image$barcode,]
#prueba_ <- prueba_ %>% mutate(barcode = gsub("_5", "", barcode))
##order barcodes
prueba_$origbarcode <- rownames(prueba_)
rownames(prueba_) <- prueba_$barcode
target <- c(rownames(M8_F2_1C_image))
prueba_ <- prueba_[match(target, prueba_$barcode),]
##
rownames(prueba_) <- prueba_$origbarcode
prueba_$barcode <- NULL
prueba_$origbarcode <- NULL

write.csv(prueba_, "./data/data/image_info/M8_coordinates.csv")



write.csv(M8_count, "data/data/image_info/M8_count.csv")
write.csv(M8, "data/data/image_info/M8.csv")








