## SCRIPT: Spatial object creation visium  2 data BM project

## 20.10.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)

#Variables---------------------------------
DIR_ROOT <- file.path(getwd())
DIR_DATA <- file.path(DIR_ROOT, "data/2.data/")

#Data--------------------------------------

samples <- dir(path = DIR_DATA)
# Create each individual Seurat object for every sample
lista <- c(samples)
names(lista) <- samples
prueba <-c()

for (i in lista){
  a <- Load10X_Spatial(data.dir = paste0(DIR_DATA, i),
                          filename = "filtered_feature_bc_matrix.h5",
                          assay = "Spatial",
                          slice = i,
                          filter.matrix = TRUE)
  prueba[[length(prueba) + 1]] <- a
}
names(prueba) <- samples

###modify image data
for (i in 1:length(prueba)) {
  #https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
  a <- prueba[[i]]
  b <- names(prueba[i])
  a@images[[b]]@coordinates[["tissue"]] <- as.integer(a@images[[b]]@coordinates[["tissue"]])
  a@images[[b]]@coordinates[["row"]] <- as.integer(a@images[[b]]@coordinates[["row"]])
  a@images[[b]]@coordinates[["col"]] <- as.integer(a@images[[b]]@coordinates[["col"]])
  a@images[[b]]@coordinates[["imagerow"]] <- as.integer(a@images[[b]]@coordinates[["imagerow"]])
  a@images[[b]]@coordinates[["imagecol"]] <- as.integer(a@images[[b]]@coordinates[["imagecol"]])
  saveRDS(a,file = paste0("./objects/sp/second/",names(lista[i]),".rds"))
}


#####################QC#############################
