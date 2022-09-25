## SCRIPT: Spatial object creation visium data BM project

## 15.09.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)

#Data--------------------------------------
M1_tib_1A <- Load10X_Spatial(data.dir = "./data/M1_tib_1A/outs/",
                      filename = "filtered_feature_bc_matrix.h5",
                      assay = "Spatial",
                      slice = "M1_tib_1A",
                      filter.matrix = TRUE)

M1_fem_2B <- Load10X_Spatial(data.dir = "./data/spatial_no_optimal_results/M1_fem_2B/outs/",
                             filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "M1_fem_2B",
                             filter.matrix = TRUE)

M3_fem_1C <- Load10X_Spatial(data.dir = "./data/spatial_no_optimal_results/M3_fem_1C/outs/",
                             filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "M3_fem_1C",
                             filter.matrix = TRUE)

M3_tib_1A <- Load10X_Spatial(data.dir = "./data/spatial_no_optimal_results/M3_tib_1A/outs/",
                             filename = "filtered_feature_bc_matrix.h5",
                             assay = "Spatial",
                             slice = "M3_tib_1A",
                             filter.matrix = TRUE)

#####add area info
######add area data
area_M1_tib_1A <- read.csv( "./data/mahtab_cloupe/M1_tibia_1A.csv")
area_M1_tib_1A <- as.vector(area_M1_tib_1A$Area)
M1_tib_1A@meta.data["area"] <- as.factor(area_M1_tib_1A)

area_M1_fem_2B <- read.csv( "./data/mahtab_cloupe/M1_femur_2B.csv")
area_M1_fem_2B <- as.vector(area_M1_fem_2B$Area)
M1_fem_2B@meta.data["area"] <- as.factor(area_M1_fem_2B)

area_M3_tib_1A <- read.csv( "./data/mahtab_cloupe/M3_tibia_1A.csv")
area_M3_tib_1A <- as.vector(area_M3_tib_1A$Area)
M3_tib_1A@meta.data["area"] <- as.factor(area_M3_tib_1A)

area_M3_fem_1C <- read.csv( "./data/mahtab_cloupe/M3_femur_1C.csv")
area_M3_fem_1C <- as.vector(area_M3_fem_1C$Area)
M3_fem_1C@meta.data["area"] <- as.factor(area_M3_fem_1C)

###########################################################

lista <- c(M1_tib_1A, M1_fem_2B, M3_fem_1C, M3_tib_1A)
names(lista) <- c("M1_tib_1A","M1_fem_2B","M3_fem_1C","M3_tib_1A")

for (i in 1:length(lista)) {
  #https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
  a <- lista[[i]]
  b <- names(lista[i])
  a@images[[b]]@coordinates[["tissue"]] <- as.integer(a@images[[b]]@coordinates[["tissue"]])
  a@images[[b]]@coordinates[["row"]] <- as.integer(a@images[[b]]@coordinates[["row"]])
  a@images[[b]]@coordinates[["col"]] <- as.integer(a@images[[b]]@coordinates[["col"]])
  a@images[[b]]@coordinates[["imagerow"]] <- as.integer(a@images[[b]]@coordinates[["imagerow"]])
  a@images[[b]]@coordinates[["imagecol"]] <- as.integer(a@images[[b]]@coordinates[["imagecol"]])
  saveRDS(a,file = paste0("./objects/sp/",names(lista[i]),".rds"))
}

##for python
a <- data.frame(M1_tib_1A@images[["M1_tib_1A"]]@coordinates)
b <- M1_tib_1A@images[["M1_tib_1A"]]@image

#save objects
write.csv(a , './PRUEBA.csv', row.names = TRUE)

