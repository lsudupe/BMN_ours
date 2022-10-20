## SCRIPT: AUCell analysis cell type BM project

## 20.10.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library("AUCell")
library("dplyr")
library("GSEABase")
library("dittoSeq")
library("dichromat")
library("escape")


#Data---------------------------------
integrated_seurat <- readRDS("./objects/sp/integrated/integrated.seurat_type.rds")
single_cell_bonemarrow <- readRDS("./objects/sc/single_cell_bonemarrow.rds")


#Analysis---------------------------------
Idents(single_cell_bonemarrow) <- single_cell_bonemarrow@meta.data[["groups"]]


##DE
groups_markers <- Seurat::FindAllMarkers(object = single_cell_bonemarrow, 
                                                  assay = "RNA",
                                                  slot = "data",
                                                  verbose = TRUE, 
                                                  only.pos = TRUE)
#Filter
markers_area <- subset(groups_markers, p_val_adj < 0.05 & 0.5 < avg_log2FC)

#Top5
markers_area <- markers_area %>%
  group_by(cluster) %>%
  top_n(n = 50,
        wt = avg_log2FC)

#create a list
cell_types <- as.list(as.vector(unique(markers_area$cluster)))
names(cell_types) <- as.vector(unique(markers_area$cluster))

#add genes into the list
for (i in 1:length(cell_types)){
  a <- cell_types[[i]]
  markers_area$cluster <- as.character(markers_area$cluster)
  b <- markers_area[markers_area$cluster %in% c(a), ]
  b <- as.vector(b$gene)
  cell_types[[i]] <- b
}

#create a list
genesets <- as.list(as.vector(unique(markers_area$cluster)))
names(genesets) <- as.vector(unique(markers_area$cluster))

#create genesets
for (i in 1:length(genesets)){
  a <- genesets[[i]]
  markers_area$cluster <- as.character(markers_area$cluster)
  b <- markers_area[markers_area$cluster %in% c(a), ]
  b <- as.vector(b$gene)
  b <- GeneSet(b, setName=paste0("geneSet",names(genesets[i])))
  genesets[[i]] <- b
}

geneSets <- GeneSetCollection(genesets)

#####AUCell
x<- integrated_seurat
x@assays[["SCT"]]@counts <- x@assays[["SCT"]]@data
matrix <- x@assays[["SCT"]]@counts
matrix <- as.matrix(matrix)
###### AUC score 
cells_rankings <- AUCell_buildRankings(matrix, nCores=1)#, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
#extract AUC values
auc_per_cell_all <- as.data.frame(t(getAUC(cells_AUC)))
##save meta
x <- AddMetaData(x, auc_per_cell_all)


#plots

set <- c("geneSetHSC_PSC", "geneSetIC", "geneSetMSC", "geneSetNC", "geneSetEC")

for (i in set){
  pdf(file.path("./results/Aucell/",filename = paste(i,"ridge.pdf",sep="")))
  print(dittoRidgePlot(x, i, group.by = "type"))
  dev.off()
}

for (i in set){
  pdf(file.path("./results/Aucell/",filename = paste(i,"ridge_2.pdf",sep="")))
  print(ridgeEnrichment(x@meta.data, gene.set = i, group = "orig.ident", facet = "area", add.rug = TRUE))
  dev.off()
}


library("GiNA")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
b <- c(0, 0.5)
for (i in set){
  p1 <- SpatialFeaturePlot(x, features = i, combine = FALSE)
  #fix.p1 <- scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98",limits = b,breaks=b, labels = c("min", "max"))
  fix.p1 <- scale_fill_gradientn(colours=myPalette(100),breaks=b, labels = c("min", "max"),limits = b)
  p2 <- lapply(p1, function (x) x + fix.p1)
  
  pdf(file.path("./results/Aucell/",filename = paste(i,"spatial.pdf",sep="")))
  print(CombinePlots(p2))
  dev.off()
}
