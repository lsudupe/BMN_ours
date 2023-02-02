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
samples <- samples[! samples %in% c("M1_tib_1A", "M3_tib_2A" ,"untitled folder","BM_human_AP-B08805")]

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
  ## subset data
  Seurat::Idents(object = a) <- a@meta.data[["area"]]
  a <- subset(x =a, idents = c("bone", "bone_marrow"))
  a@meta.data[["area"]] <- a@active.ident
  ## fix scalefactor
  #https://github.com/satijalab/seurat/issues/5614
  a@images[[i]]@scale.factors$lowres = a@images[[i]]@scale.factors$hires
  ## fix image
  #https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
  a@images[[i]]@coordinates[["tissue"]] <- as.integer(a@images[[i]]@coordinates[["tissue"]])
  a@images[[i]]@coordinates[["row"]] <- as.integer(a@images[[i]]@coordinates[["row"]])
  a@images[[i]]@coordinates[["col"]] <- as.integer(a@images[[i]]@coordinates[["col"]])
  a@images[[i]]@coordinates[["imagerow"]] <- as.integer(a@images[[i]]@coordinates[["imagerow"]])
  a@images[[i]]@coordinates[["imagecol"]] <- as.integer(a@images[[i]]@coordinates[["imagecol"]])
  ######add object to list
  prueba[[length(prueba) + 1]] <- a
  
}
names(prueba) <- samples

## separate list
list2env(prueba,envir=.GlobalEnv)

## add sample info
M1_fem_1C@meta.data[["type"]] <- "MM"
M2_F_2B@meta.data[["type"]] <- "MM"
M3_F_1C@meta.data[["type"]] <- "control"
M3_fem_1C@meta.data[["type"]] <- "control"
M8_F2_1C@meta.data[["type"]] <- "MM"
M9_F2_1C@meta.data[["type"]] <- "MM"

## combine
combined <- merge(M1_fem_1C, y = c(M2_F_2B, M3_F_1C, M3_fem_1C, M8_F2_1C, M9_F2_1C ), 
                  add.cell.ids = c("M1_fem_1C", "M2_F_2B", "M3_F_1C", "M3_fem_1C", "M8_F2_1C", "M9_F2_1C"), project = "BM")

saveRDS(combined, "./objects/sp/combined.rds")
###########################################################

##for python
a <- data.frame(M1_tib_1A@images[["M1_tib_1A"]]@coordinates)
b <- M1_tib_1A@images[["M1_tib_1A"]]@image

#save objects
write.csv(a , './PRUEBA.csv', row.names = TRUE)

