## SCRIPT: AUCell analysis cell type BM project

## 20.10.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library("AUCell")
library("dplyr")
library("GSEABase")


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


###### Create genesets
nine.sig <- GeneSet(nine.no.b, setName="geneSet9")