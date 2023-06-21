## SCRIPT: Plotting script Bone Marrow project

## 21.06.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(STutility)

# Data
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

b <- SetIdent(se, value = se@meta.data[["clustering"]])
b <- SubsetSTData(b, idents = c("5", "6","7"))
b@meta.data[["hotspots"]] <- b@active.ident

b <- SetIdent(b, value = b@meta.data[["name"]])
subset <- SubsetSTData(b, idents = c("M1_fem_1C", "M2_F_2B", "M8_F2_1C", "M9_F2_1C"))

#spatial
##divide by sample
Idents(object = subset) <- "name"
name <- unique(subset@meta.data[["name"]])

objects <- c()
subset <- SetIdent(subset, value = subset@meta.data[["name"]])

for (i in name){
  a <- SubsetSTData(subset, idents = i)
  ## add object to list
  objects[[length(objects) + 1]] <- a
}

names(objects) <- name
## separate list
list2env(objects,envir=.GlobalEnv)

###create a loop, 
###separate each object in 3 hotspots
###plot each variable of interest in each object
###enrichment scores

for (i in 1:length(objects)){
  a <- objects[[i]]
  a <- SetIdent(a, value = a@meta.data[["hotspots"]])
  b <- SubsetSTData(b, idents = c("5", "6","7"))
  
  
}






