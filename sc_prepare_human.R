## SCRIPT: Singel cell reference data HUMAN preparing BM project

## 14.06.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2) 
library(tidyverse)
library(base)
library(harmony)


#Data--------------------------------------
# Read the raw counts data
raw_counts <- read.table("./data/single-cell/human/neutrophil/GSE117131_Zhu_SC-10X_human_raw_counts.txt", header = TRUE, row.names = 1)

# Read the annotation data
annotations <- read.table("./data/single-cell/human/neutrophil/GSE117131_Zhu_SC-10X_human_annotation.txt", header = TRUE, row.names = 1)

##check
all(colnames(raw_counts) %in% rownames(annotations)) # should return TRUE
# Create the Seurat object with the raw counts and the metadata
neutro <- CreateSeuratObject(counts = raw_counts, project = "Zhu_SC-10X", meta.data = annotations)
# Normalize the data
seurat_object <- NormalizeData(neutro)
# Find the most highly expressed genes
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# Perform a scaling of the data
seurat_object <- ScaleData(seurat_object)
# Determine QC metrics
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
# Visualize QC metrics
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
saveRDS(seurat_object, "./objects/sc/human/neutro.rds")

##pc
pc_mm <- readRDS("./data/single-cell/human/pc_mm/MM_multiome_list.rds")

m1 <- pc_mm[["MM1"]]
# Set the default assay to RNA
counts <- m1@assays[["RNA"]]@counts
meta <- m1@meta.data
m1 <- CreateSeuratObject(counts = counts,meta.data = meta)
Seurat::Idents(object = m1) <- m1@meta.data[["orig.ident"]]

neutro <- readRDS("./objects/sc/human/neutro.rds")
all <- readRDS("./data/single-cell/human/bone_marrow/disco_bone_marrow_v01 (1).rds")

##transform all data
counts <- all@assays[["RNA"]]@counts
meta <- all@meta.data
ref <- CreateSeuratObject(counts = counts,meta.data = meta)
##transform all data fin

#####combine data
neutro@meta.data[["orig.ident"]] <- "neutro"
Seurat::Idents(object = neutro) <- neutro@meta.data[["orig.ident"]]
ref@meta.data[["orig.ident"]] <- "ref"
Seurat::Idents(object = ref) <- ref@meta.data[["ct"]]

combined_sc <- merge(ref, y = c(neutro, m1), add.spot.ids = c("ref", "neutro", "m1"), project = "BM")

combined_sc@meta.data[["split"]] <- combined_sc@active.ident
Seurat::Idents(object = combined_sc) <- combined_sc@meta.data[["orig.ident"]]
harmony <- combined_sc

meta <- harmony@meta.data
meta <- meta[1:18]
harmony@meta.data <- meta
saveRDS(harmony, "./objects/sc/integrated/integrated_sc_human_merge.rds")
harmony <- readRDS("./objects/sc/integrated/integrated_sc_human_merge.rds")

## Harmony
harmony_i <- harmony %>%
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(npcs = 15) %>%
  RunHarmony(assay.use="RNA",reduction = "pca", dims = 1:15, group.by.vars = "orig.ident") %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution=0.7) %>%
  RunUMAP(reduction = "harmony", dims = 1:15, n.epochs = 1e3) 

saveRDS(harmony_i, "./objects/sc/integrated/integrated_sc_human_harmony.rds")
saveRDS(subset, "./objects/sc/integrated/single_cell_bonemarrow_human_all_groups_harmony.rds")

###extract neutrophile marker genes
group_of_interest <- "neutro"

# Find markers for that specific group compared to all other cells
markers <- FindMarkers(harmony, ident.1 = group_of_interest)

# Extract the top 100 genes (sorted by adjusted p-value)
top_genes <- head(markers, 50)$gene

# Print the top genes
print(top_genes)



