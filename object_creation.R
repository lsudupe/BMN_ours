## SCRIPT: Spatial object creation visium data BM project

## 15.09.22 Laura Sudupe , git @lsudupe

## Libraries
library(Seurat)
library(tidyverse)



## Read data
#Variables---------------------------------
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "data/data/")

## samples
samples <- dir(path = DIR_DATA)
samples <- samples[! samples %in% c("M1_tib_1A", "M3_tib_2A" )]

## Create each individual Seurat object for every sample
lista <- c(samples)
names(lista) <- samples
prueba <-c()

for (i in lista){
  a <- Load10X_Spatial(data.dir = paste0(DIR_DATA, i),
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial",
                       slice = i,
                       filter.matrix = TRUE)
  
  ## add area data
  area_a <- read.csv(paste0(DIR_DATA, i, "/area.csv"))
  area_a <- as.vector(area_a$area)
  a@meta.data["area"] <- as.factor(area_a)
  a@meta.data["orig.ident"] <- i
  ######subset data
  Seurat::Idents(object = a) <- a@meta.data[["area"]]
  a <- subset(x =a, idents = c("bone", "bone_marrow"))
  a@meta.data[["area"]] <- a@active.ident
  ######add object to list
  prueba[[length(prueba) + 1]] <- a
  
}
names(prueba) <- samples

## separate list
list2env(prueba,envir=.GlobalEnv)

## combine
combined <- merge(M1_fem_1C, y = c(M2_F_2B, M3_F_1C, M3_fem_1C, M8_F2_1C, M9_F2_1C ), 
                  add.cell.ids = c("M1_fem_1C", "M2_F_2B", "M3_F_1C", "M3_fem_1C", "M8_F2_1C", "M9_F2_1C"), project = "BM")

saveRDS(combined, "./objects/sp/combined.rds")
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

