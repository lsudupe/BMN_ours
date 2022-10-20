## SCRIPT: Mahtab division cell type endogram BM project

## 19.10.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
#BiocManager::install("HGC")
library(HGC)
#https://www.datanovia.com/en/lessons/examples-of-dendrograms-visualization/
library(factoextra)
#Data---------------------------------
integrated_seurat <- readRDS("./objects/sp/integrated/integrated.seurat_type.rds")


#Analysis---------------------------------


####find tree seurat
Idents(integrated_seurat) <- integrated_seurat@meta.data[["area"]]
integrated_seurat <- FindClusteringTree(integrated_seurat, graph.type = "KNN")
seurat_tree <- integrated_seurat@graphs$ClusteringTree
seurat_data <- as.data.frame(integrated_seurat@meta.data)
seurat.labels <- data.frame(area = seurat_data$area,
                            prueba = seurat_data$orig.ident)

seurat_tree$height = log(seurat_tree$height + 1)
seurat_tree$height = log(seurat_tree$height + 1)

seurat_tree[["labels"]] <- integrated_harmony@meta.data[["area"]]

pdf(file.path("./results/endogram/",filename = "endo.pdf"))
print(HGC.PlotDendrogram(tree = seurat_tree,
                         k = 5, plot.label = TRUE,
                         labels = seurat.labels))
dev.off()

pdf(file.path("./results/endogram/",filename = "seurat_spatial.pdf"))
print(SpatialDimPlot(integrated_seurat, label = T, combine = F))
dev.off()

a <- HGC.PlotDendrogram(tree = seurat_tree,
                        k = 5, plot.label = TRUE,
                        labels = seurat.labels)
library(plotly)
library(ggplot2)
library(ggdendro)

p <- ggdendrogram(seurat_tree, rotate = FALSE, size = 2)
print(ggplotly(p)) 


library(BuenColors)
nuria <- jdb_palette("brewer_spectra")
library(RColorBrewer)
a <- rev(brewer.pal(n = 19, name = "RdBu"))
b <- c(nuria + a)
lista <- c(nuria, a)

res <- hcut(USArrests, k = 4, stand = TRUE)

pdf(file.path("./results/endogram/",filename = "seurat.pdf"))
print(fviz_dend(seurat_tree, rect = TRUE, k_colors = lista))
dev.off()


data("USArrests")
df <- scale(USArrests)
res <- hcut(USArrests, k = 4, stand = TRUE)
