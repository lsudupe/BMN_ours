

library(Seurat)


meta <- read.table("./data/single-cell/human/GSE117131_Zhu_SC-10X_human_annotation.txt")
result <- meta[-1]
row.names(result) <- meta$Barcode


count <- read.table("./data/single-cell/human/GSE117131_Zhu_SC-10X_human_raw-counts.txt")

object <- CreateSeuratObject(
  prueba2,
  project = "bm_human",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = result
)


object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)

object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(object)
object <- ScaleData(object, features = all.genes)
object <- RunPCA(object, features = VariableFeatures(object = object))

object <- FindNeighbors(object, dims = 1:10)
object <- FindClusters(object, resolution = 0.5)
object <- RunUMAP(object, dims = 1:10)

DimPlot(object, reduction = "umap")

FeaturePlot(object, features = c("IFI6"))

