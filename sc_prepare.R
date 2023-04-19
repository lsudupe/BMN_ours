## SCRIPT: Singel cell reference data preparing BM project

## 23.11.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(scCustomize)
library(dplyr)
library(ggplot2) 
library('magrittr')
library(tidyverse)
library(base)
library(harmony)
library(scclusteval)

#Data--------------------------------------
load(file='./data/single-cell/RNAMagnetDataBundle/NicheMarkers10x.rda')
load(file='./data/single-cell/RNAMagnetDataBundle/NicheData10x.rda')
load(file="./data/single-cell/RNAMagnetDataBundle/NicheDataLCM.rda")
load(file="./data/single-cell/RNAMagnetDataBundle/NicheMetaDataLCM.rda")
single_cell_bonemarrow = UpdateSeuratObject(object = NicheData10x)

single_cell_bonemarrow@meta.data[["cell"]] <- single_cell_bonemarrow@active.ident
#####SUBSET IN BONE, BONE MARROW

a <- as.factor(single_cell_bonemarrow@meta.data[["cell"]])
levels(a) <- list(PSC  = "Ery/Mk prog.", 
                  PSC = "Neutro prog.",
                  PSC = "Mono prog.",
                  PSC = "Gran/Mono prog.",
                  PSC = "LMPPs",
                  Bcell = "large pre-B.",
                  PSC = "Mk prog.",
                  Erythroblasts = "Erythroblasts",
                  PSC = "Eo/Baso prog.",
                  Monocytes = "Monocytes",
                  PSC = "Ery prog.",
                  Bcell = "pro-B",
                  Tcell = "T cells",
                  Neutrophils = "Neutrophils",
                  MSC = "Adipo-CAR",
                  MSC = "Ng2+ MSCs",
                  NC = "Schwann cells",
                  MSC = "Osteoblasts",
                  Fibro = "Arteriolar fibro.",
                  EC = "Sinusoidal ECs",
                  MSC = "Osteo-CAR",
                  Bcell = "small pre-B.",
                  MSC = "Chondrocytes",
                  Fibro = "Endosteal fibro.",
                  Fibro = "Fibro/Chondro p.",
                  Fibro = "Stromal fibro.",
                  EC = "Arteriolar ECs",
                  Fibro = "Myofibroblasts",
                  MSC = "Smooth muscle",
                  DC = "Dendritic cells",
                  NK = "NK cells",
                  Bcell = "B cell"
)

#SMC, NC, Fibro_Chondro_p take out, 
single_cell_bonemarrow@meta.data[["ident"]] <- a
Seurat::Idents(object = single_cell_bonemarrow) <- single_cell_bonemarrow@meta.data[["ident"]]

single_cell_bonemarrow <- subset(x = single_cell_bonemarrow, idents = c("NC", "PSC", "Fibro"), invert = TRUE)
single_cell_bonemarrow@meta.data[["ident"]] <- single_cell_bonemarrow@active.ident


saveRDS(single_cell_bonemarrow, "./objects/heterogeneity/single_cell_bonemarrow.rds")
saveRDS(single_cell_bonemarrow, "./objects/heterogeneity/single_cell_bonemarrow_all_groups.rds")
single_cell_bonemarrow <- readRDS("./objects/heterogeneity/single_cell_bonemarrow.rds")

pdf(file.path("./results/clusters/single_cell/",filename = "subgroups_sc.pdf"))
#print(DimPlot(single_cell_bonemarrow, reduction = "tsne", label = TRUE))
print(DimPlot(single_cell_bonemarrow, group.by = c("ident"), label = TRUE, repel=TRUE))
dev.off()



############READ AZARI data
PC_MM <- readRDS("./data/single-cell/PC/scRNA_MM_PC.rds")
PC_MM <- subset(x = PC_MM, idents = c("MM_MIC"))
PC_MM@meta.data[["ident"]] <- PC_MM@active.ident


#####integrate data
##Merge them
single_cell_bonemarrow@meta.data[["orig.ident"]] <- "ref"
combined_sc <- merge(single_cell_bonemarrow, y = c( PC_MM), 
                  add.cell.ids = c("single_cell_bonemarrow", "PC_MM"), project = "BM")
combined_sc@meta.data[["split"]] <- combined_sc@active.ident
Seurat::Idents(object = combined_sc) <- combined_sc@meta.data[["orig.ident"]]
x <- combined_sc
###separate the data
seurat <- x
harmony <- x

## SEURAT
# Select the most variable features to use for integration
options(future.globals.maxSize = 2000 * 1024^2)  # set allowed size to 2K MiB
list <- SplitObject(seurat, split.by = "orig.ident")
list <- lapply(X = list, FUN = SCTransform, assay="RNA")
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
seurat_i <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

seurat_i <- RunPCA(seurat_i, verbose=FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution=0.7) 

saveRDS(seurat_i, "./objects/sc/integrated/integrated_sc_seurat.rds")

## Harmony
harmony_i <- harmony %>%
    NormalizeData() %>%
    ScaleData() %>%
    FindVariableFeatures() %>%
    RunPCA(npcs = 20) %>%
    RunHarmony(assay.use="RNA",reduction = "pca", dims = 1:20, group.by.vars = "orig.ident") %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution=0.7) %>%
    RunUMAP(reduction = "harmony", dims = 1:20, n.epochs = 1e3) 

saveRDS(harmony_i, "./objects/sc/integrated/integrated_sc_harmony.rds")
saveRDS(harmony_i, "./objects/sc/integrated/single_cell_bonemarrow_all_groups_harmony.rds")

harmony <- readRDS("./objects/sc/integrated/integrated_sc_harmony.rds")

Seurat::Idents(object = harmony) <- harmony@meta.data[["ident"]]

pdf(file.path("./results/clusters/single_cell/",filename = "groups_sc.pdf"))
#print(DimPlot(single_cell_bonemarrow, reduction = "tsne", label = TRUE))
print(DimPlot(harmony, group.by = c("ident"), label = TRUE, repel=TRUE))
dev.off()

## plots
samples <- c(seurat_i, harmony_i)
names(samples) <- c("seurat_i", "harmony_i")

for (i in 1:length(samples)){
  a <- samples[[i]]
  #umap separate in cell type
  pdf(file.path("./results/clusters/single_cell/integration/",filename = paste0("umap_orig.ident_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("orig.ident"), label = T) + ggtitle("cell type"))
  dev.off()
  #umap separate in samples
  pdf(file.path("./results/clusters/single_cell/integration/",filename = paste0("umap_split_",names(samples[i]),".pdf")))
  print(DimPlot(a, group.by = c("split"), label = T) + ggtitle("sample"))
  dev.off()
}

######Jackard index harmony clusters vs seurat
Seurat::Idents(object = seurat_i) <- seurat_i@meta.data[["split"]]
Seurat::Idents(object = harmony_i) <- harmony_i@meta.data[["split"]]

pdf(file.path("./results/jackard/",filename = "jackard_type_harmonyvsseurat_singlecell.pdf"))                                                         
PairWiseJaccardSetsHeatmap(seurat_i@active.ident,
                           harmony_i@active.ident,
                           best_match = TRUE)
dev.off()

####marker genes for MM
# Determine differentiating markers for PC_MM
Seurat::Idents(object = PC_MM) <- PC_MM@meta.data[["orig.ident"]]
PC_MM <- PrepSCTFindMarkers(PC_MM, assay = "SCT", verbose = TRUE)
MM_genes <- Seurat::FindAllMarkers(object = PC_MM, 
                                              assay = "SCT",
                                              verbose = TRUE, 
                                              only.pos = TRUE)
cluster_MIC <- subset(MM_genes, p_val_adj < 0.05 & 0.05 < avg_log2FC)
genes <- cluster_MIC[grepl("MM_MIC", cluster_MIC[,6]),]
genes <- cluster_MIC$gene

#Top5
top50 <- cluster_6_f %>%
  top_n(n = 50,
        wt = avg_log2FC)

write.csv(genes, "./data/single-cell/PC/MM_MIC_genes.csv", row.names = FALSE)

#########CELL SCORE CHECK########
PC_MM_S <- readRDS("./data/single-cell/PC/scRNA_MM_PC.rds")
DefaultAssay(PC_MM_S) <- "RNA"

PC_MM_S <- PC_MM_S %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(npcs = 20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution=0.7) %>%
  RunUMAP(reduction = "pca", dims = 1:20, n.epochs = 1e3) 


###cell score plots##
PC_MM@meta.data[["Phase"]]
DimPlot(PC_MM_S, group.by = c("Phase"), label = TRUE, repel=TRUE)

####marker genes for MM
# Determine differentiating markers for PC_MM
Seurat::Idents(object = PC_MM_S) <- PC_MM_S@meta.data[["Phase"]]
#PC_MM_S <- PrepSCTFindMarkers(PC_MM_S, assay = "SCT", verbose = TRUE)
MM_genes <- Seurat::FindAllMarkers(object = PC_MM_S, 
                                   assay = "RNA",
                                   verbose = TRUE, 
                                   only.pos = TRUE)
cluster <- subset(MM_genes, p_val_adj < 0.05 & 0.05 < avg_log2FC)
genes_G2M <- cluster[grepl("G2M", cluster[,6]),]


genes_G1 <- cluster[grepl("G1", cluster[,6]),]
genes_G1 <- cluster$gene

genes_S <- cluster[grepl("S", cluster[,6]),]
genes_S <- cluster$gene

#Top5
top50_G2M <- genes_G2M %>%
  top_n(n = 50,
        wt = avg_log2FC)
genes_G2M <- top50_G2M$gene

top50_G1 <- genes_G1 %>%
  top_n(n = 50,
        wt = avg_log2FC)
genes_G1 <- top50_G1$gene

top50_S <- genes_S %>%
  top_n(n = 50,
        wt = avg_log2FC)
genes_S <- top50_S$gene

intersect <- intersect(genes_G2M,genes_G1)

###Add module score
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")
se <- AddModuleScore(se,
                       features = list(genes_G1),
                       name="G1")

se <- AddModuleScore(se,
                     features = list(genes_G2M),
                     name="G2M")

se <- AddModuleScore(se,
                     features = list(genes_G2M),
                     name="S")

library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(11,"Spectral")
color <- rev(color)

FeatureOverlay(se, features = c("G11"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = color)
FeatureOverlay(se, features = c("G2M1"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
               value.scale = "all" ,cols = color)
FeatureOverlay(se, features = c("S1"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
               value.scale = "all" ,cols = color)

###Add module score FIN





