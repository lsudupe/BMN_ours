## SCRIPT: Heterogeneity BONE MARROW

## 05.12.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------

library(Seurat)
library(Matrix)
library(genefilter)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
color <- rev(brewer.pal(11,"Spectral"))

## Data
M1 <- readRDS("./objects/sp/st_s1_module.rds")
M2 <- readRDS("./objects/sp/st_s2_module.rds")
M8 <- readRDS("./objects/sp/st_s8_module.rds")
M9 <- readRDS("./objects/sp/st_s9_module.rds")

identifiers <- c("M1", "M2", "M8", "M9")
# Use MergeSTData to merge the objects
se.merged <- MergeSTData(
  x = M1,
  y = list(M2, M8, M9),
  add.spot.ids = identifiers,
  merge.data = TRUE
)

# Find immunoglobulin genes by pattern
immunoglobulin_genes <- grep("^Igh", rownames(se.merged), value = TRUE)
# Remove these genes from the matrix
se.merged <- se.merged[!rownames(se.merged) %in% immunoglobulin_genes, ]

metaData <- se.merged@meta.data
metaData <- data.frame(row.names = rownames(metaData),
                       sampleID = metaData$name,
                       modules = metaData$community)

st <- SCTransform(se.merged, assay = "RNA", new.assay.name = "SCT_novars")

st <- SCTransform(se.merged, assay = "RNA", new.assay.name = "SCT_novars_sample", vars.to.regress = "name")
#st <- SCTransform(se.merged, assay = "RNA", new.assay.name = "SCT_novars_pc", vars.to.regress = "name")

# Step1: aggregate counts by community
Seurat::Idents(object = st) <- st@meta.data[["community"]]
cts <- AggregateExpression(st, 
                           group.by = c("community"),
                           assays = 'SCT_novars',
                           #assays = 'SCT_novars_sample',
                           slot = "data",
                           return.seurat = FALSE)

aggregated_counts <- cts$SCT_novars
#aggregated_counts <- cts$SCT_novars_sample
head(aggregated_counts)
row_totals <- rowSums(aggregated_counts)
filtered_data <- aggregated_counts[row_totals >= 10, ]
row_totals <- rowSums(filtered_data)
quantile_99 <- quantile(filtered_data, probs = 0.99)
#filtered_data <- filtered_data[row_totals <= 500, ]
head(filtered_data)

final_mat <- NormalizeData(filtered_data)
final_mat <- ScaleData(final_mat)
head(final_mat)

final_mat <- NormalizeData(filtered_data)
final_mat <- ScaleData(final_mat)

## top var genes
### standar error
se <- function(x) sd(x)/sqrt(length(x))
t_final_mat <- Matrix::t(final_mat)
se_mat <- apply(t_final_mat, 2, se)
se_mat <- se_mat[order(se_mat, decreasing = T)]
se_mat_top100 <- names(sort(se_mat, decreasing = TRUE))[1:100]

## pca
pca<- prcomp(t(expression_mat_sub),center = TRUE, scale. = TRUE) 
# check the order of the samples are the same.
all.equal(rownames(pca$x), colnames(expression_mat_sub))
# extract pca 1 and pca2
PC1_and_PC2<- data.frame(PC1=pca$x[,1], PC2= pca$x[,2])
PC1_and_PC2$sample <- c("M1","M1", "M2", "M2","M2","M2","M2", "M8","M8","M9","M9","M9")
head(PC1_and_PC2)

ggplot(PC1_and_PC2, aes(x=PC1, y=PC2)) + 
geom_point(aes(color = sample)) +
theme_bw(base_size = 14) 

ggplot(PC1_and_PC2, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color = rownames(PC1_and_PC2))) +
  theme_bw(base_size = 14) 

## data
#scale.data <- st@assays$SCT_novars_sample@scale.data #826 spots, 3000 genes
#dim(scale.data)
scale.data <- st@assays$SCT_novars@scale.data #826 spots, 3000 genes
dim(scale.data)
scale.data_se <- scale.data[rownames(scale.data) %in% names(se_mat)[1:100],] #826 spots 76 genes
dim(scale.data_se)


metaData <- st@meta.data
metaData <- data.frame(row.names = rownames(metaData),
                       modules = metaData$community)


## heatmap
library(ComplexHeatmap)
Heatmap(scale.data_se, show_column_names = F, row_names_gp = gpar(fontsize=6),cluster_columns = FALSE,
        top_annotation = HeatmapAnnotation(df = metaData,
                                           col = list(modules = c("S1_M1" = "#1f77b4",  # blue
                                                                  "S1_M2" = "#e377c2",  
                                                                  "S2_M1" = "#ff7f0e",  # orange
                                                                  "S2_M2" = "#2ca02c",  # green
                                                                  "S2_M3" = "#d62728",  # red
                                                                  "S2_M4" = "#9467bd",  # purple
                                                                  "S2_M5" = "#731e59",
                                                                  "S8_M1" = "#8c564b",  # brown
                                                                  "S8_M2" = "#e377c2",  # pink
                                                                  "S9_M1" = "#7f7f7f",  # grey
                                                                  "S9_M2" = "#bcbd22",  # yellow-green
                                                                  "S9_M3" = "#17becf"))))   # cyan))))


module_colors <- c(S1_M1 = "#1f77b4",  # blue
                   S1_M2 = "#e377c2",  # pink
                   S2_M1 = "#ff7f0e",  # orange
                   S2_M2 = "#2ca02c",  # green
                   S2_M3 = "#d62728",  # red
                   S2_M4 = "#9467bd",  # purple
                   S2_M5 = "#731e59",  # dark purple
                   S8_M1 = "#8c564b",  # brown
                   S8_M2 = "#e377c2",  # pink
                   S9_M1 = "#7f7f7f",  # grey
                   S9_M2 = "#bcbd22",  # yellow-green
                   S9_M3 = "#17becf")  # cyan

# Use the module_colors named vector directly
top_annotation <- HeatmapAnnotation(modules = module_colors)

se_mat_top500 <- names(sort(se_mat, decreasing = TRUE))[1:100]

final_mat_se <- final_mat[se_mat_top500, ]
final_mat_se <- final_mat[se_mat_top100, ] #pseudobulk ordered most var genes
dim(final_mat_se)

prueba <- final_mat_se
prueba <- NormalizeData(prueba)
prueba <- ScaleData(prueba)
# Draw the heatmap
suppressMessages(
  Heatmap(prueba, show_column_names = TRUE, 
          row_names_gp = gpar(fontsize=6), cluster_columns = FALSE,
          top_annotation = top_annotation)
)

#noRegressGenes <- se_mat_top100
nameRegressGenes <- se_mat_top100
pcRegreesGenes <- se_mat_top100

top5 <- names(se_mat)[1:20]
top5

x <- "Igkc"


m1 <- FeatureOverlay(M1, features = c(x),pt.size = 1.3, col=color)
m2 <- FeatureOverlay(M2, features = c(x),pt.size = 1.3, col=color)
m3 <- FeatureOverlay(M8, features = c(x),pt.size = 1.3, col=color)
m4 <- FeatureOverlay(M9, features = c(x),pt.size = 1.3, col=color)

grid.arrange(m1, m2, m3, m4, ncol = 2)


VlnPlot(st, features = x, group.by = "community", assay = "SCT")



