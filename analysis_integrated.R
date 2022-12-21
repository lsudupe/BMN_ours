## SCRIPT: Integration spatial data BMN project

## 20.09.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2) 
library(harmony)
library('magrittr')
library(tidyverse)
library(base)
library(scclusteval)

#Data--------------------------------------
combined  <- readRDS("./objects/sp/second/combined_filtered.rds")
x <- combined

x.image <- x@images
x@images[["M1_tib_1A"]]<- NULL
x@images[["M1_fem_1C"]]<- NULL
x@images[["M3_tib_2A"]]<- NULL
x@images[["M3_fem_1C"]]<- NULL

###separate the data
seurat_resolution <- x
harmony_resolution <- x

######seurat#####
list <- SplitObject(seurat_resolution, split.by = "type")
list <- lapply(X = list, FUN = SCTransform, assay="Spatial")
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
seurat_resolution <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

seurat_resolution <- RunPCA(seurat_resolution,npcs = 15, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:15, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:15) %>%
  FindClusters(resolution=0.7) 

seurat_resolution@images <- x.image

#a <- as.factor(integrated_seurat@meta.data[["area"]])
#levels(a) <- list(BM_SOC  = "bone marrow (SOC)", 
 #                 AC_mix = "AC mix",
  #                DT = "DT",
   #               cortical_bone = "cortical bone",
    #              BM_B = "bone marrow/bone",
     #             NS = "NS",
      #            adipocyte_BM = "adipocyte in bone marrow",
       #           adipocyte_M = "adipocyte in muscle",
#                  AC = "AC",
 #                 muscle = "muscle",
  #                endosteal_B = "endosteal bone",
   #               periosteal_B = "periosteal bone",
    #              trabecular_bone_in_BM = "trabecular bone in bone marrow",
     #             SOC_GP = "SOC-GP",
      #            GP_POC = "GP-POC",
       #           BM = "bone marrow",
        #          B = "bone")
#integrated_seurat@meta.data[["area"]] <- a
saveRDS(seurat_resolution, "./objects/sp/integrated/second/integrated.seurat_type.rds")

######harmony#####
harmony_resolution <- harmony_resolution %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(npcs = 15) %>%
  RunHarmony(assay.use="Spatial",reduction = "pca", dims = 1:15, group.by.vars = "type") %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution=0.7) %>%
  RunUMAP(reduction = "harmony", dims = 1:15, n.epochs = 1e3) 

harmony_resolution@images <- x.image
saveRDS(harmony_resolution, "./objects/sp/integrated/second/integrated.harmony_type.rds")

integrated_seurat <- readRDS("./objects/sp/integrated/second/integrated.seurat_type.rds")
integrated_harmony <- readRDS("./objects/sp/integrated/second/integrated.harmony_type.rds")

######Jackard index harmony clusters vs seurat
pdf(file.path("./results/jackard/second",filename = "jackard_type.pdf"))                                                         
PairWiseJaccardSetsHeatmap(integrated_seurat@active.ident,
                           integrated_harmony@active.ident)
dev.off()

######Jackard index seurat clusters vs area
clusters <- integrated_seurat
clusters <- SetIdent(clusters, value = clusters@meta.data[["seurat_clusters"]])
area <- integrated_seurat
area <- SetIdent(area, value = area@meta.data[["area"]])

pdf(file.path("./results/jackard",filename = "jackard_seurat_clusters_vs_area.pdf"))                                                         
PairWiseJaccardSetsHeatmap(clusters@active.ident,
                           area@active.ident)
dev.off()

######Jackard index harmony clusters vs area
clusters <- integrated_harmony
clusters <- SetIdent(clusters, value = clusters@meta.data[["seurat_clusters"]])
area <- integrated_seurat
area <- SetIdent(area, value = area@meta.data[["area"]])

pdf(file.path("./results/jackard",filename = "jackard_harmony_clusters_vs_area.pdf"))                                                         
PairWiseJaccardSetsHeatmap(clusters@active.ident,
                           area@active.ident)
dev.off()

samples <- c(integrated_seurat, integrated_harmony)
names(samples) <- c("integrated_seurat", "integrated_harmony")

for (i in 1:length(samples)){
  a <- samples[[i]]
  #umap separate in cell type
  pdf(file.path("./results/clusters/integration/second/",filename = paste0("umap_clusters_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("seurat_clusters"), label = T) + ggtitle("cell type"))
  dev.off()
  #umap separate in samples
  pdf(file.path("./results/clusters/integration/second/",filename = paste0("umap_samples_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("orig.ident"), label = T) + ggtitle("sample"))
  dev.off()
  #umap separate in samples
  #pdf(file.path("./results/clusters/integration/",filename = paste0("umap_area_",names(samples[i]),".pdf")))
  #print(DimPlot(a, group.by = c("area"),repel=TRUE, label = T) + ggtitle("sample"))
  #dev.off()
  #umap separate in samples
  pdf(file.path("./results/clusters/integration/second/",filename = paste0("umap_type_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("type"), label = T) + ggtitle("sample"))
  dev.off()
  a <- SetIdent(a, value = a@meta.data[["seurat_clusters"]])
  #spatial umap 
  pdf(file.path("./results/clusters/integration/second/",filename = paste0("cluster_spatial_",names(samples[i]),".pdf")))
  print(SpatialDimPlot(a, combine = FALSE,label.size = 1.5, label = T, crop = TRUE, pt.size.factor = 2.5))
  dev.off()
  #a <- SetIdent(a, value = a@meta.data[["area"]])
  #spatial umap 
  #pdf(file.path("./results/clusters/integration/",filename = paste0("area_spatial_",names(samples[i]),".pdf")))
  #print(SpatialDimPlot(a, combine = FALSE,label.size = 1.5, label = T))
  #dev.off()
}


###PLOTS###
######quality metrics per area
pdf(file.path("./results/clusters/integration/",filename = "Number of spot per area.pdf"))
integrated_seurat@meta.data%>% 
  ggplot(aes(x=area, fill=orig.ident)) + 
  geom_bar(alpha=0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NSpots")
dev.off()

pdf(file.path("./results/clusters/integration/",filename = "Number of spot per area tissue.pdf"))
integrated_seurat@meta.data%>% 
  ggplot(aes(x=orig.ident, fill=area)) + 
  geom_bar(alpha=0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NSpots")
dev.off()

###Markers

x <- integrated_seurat
x <- SetIdent(x, value = x@meta.data[["area"]])
markers_seurat_area <- Seurat::FindAllMarkers(object = x, 
                                    assay = "Spatial",
                                    slot = "data",
                                    verbose = TRUE, 
                                    only.pos = TRUE)

x <- SetIdent(x, value = area@meta.data[["seurat_clusters"]])
markers_seurat_clusters <- Seurat::FindAllMarkers(object = x, 
                                         assay = "Spatial",
                                         slot = "data",
                                         verbose = TRUE, 
                                         only.pos = TRUE)

#saveRDS(markers_x, "./results/marker_genes/M1_tib_1A.rds")

#Filter
markers_area <- subset(markers_seurat_area, p_val_adj < 0.05 & 0.5 < avg_log2FC)
markers_seurat <- subset(markers_seurat_clusters, p_val_adj < 0.05 & 0.5 < avg_log2FC)

#Top5
top20_area <- markers_area %>%
  group_by(cluster) %>%
  top_n(n = 20,
        wt = avg_log2FC)

top20_seurat <- markers_seurat %>%
  group_by(cluster) %>%
  top_n(n = 20,
        wt = avg_log2FC)

write.csv(top20_area, "./results/DE/emma/top20_genes_area.csv", row.names=FALSE)


cell_types <- as.list(as.vector(unique(top20_area$cluster)))
names(cell_types) <- as.vector(unique(top20_area$cluster))

for (i in 1:length(cell_types)){
  a <- cell_types[[i]]
  top20_area$cluster <- as.character(top20_area$cluster)
  b <- top20_area[top20_area$cluster %in% c(a), ]
  b <- as.vector(b$gene)
  cell_types[[i]] <- b
}

integrated_seurat <- readRDS("./objects/sp/integrated/integrated.seurat_type.rds")


####AddmoduleScore with Mahtab annotation
DefaultAssay(integrated_seurat) <- "SCT"
for (i in 1:length(cell_types)){
  a <- cell_types[[i]]
  c <- paste0(names(cell_types[i]), "1")
  integrated_seurat <- AddModuleScore(integrated_seurat, features = cell_types[i], name = names(cell_types[i]))
  #spatial 
  pdf(file.path("./results/marker_genes/areas_Mahtab/seurat/",filename = paste0("spatial_seurat_",c,".pdf")))
  print(SpatialFeaturePlot(integrated_seurat,features=c,combine = FALSE))
  dev.off()
}

DoHeatmap(integrated_seurat,assay = "SCT", group.by = integrated_seurat@meta.data[["type"]],
          features = cell_types[[1]])
DoHeatmap(integrated_seurat,assay = "SCT", features = top20_seurat$gene)
