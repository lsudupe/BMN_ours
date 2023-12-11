## SCRIPT: Heterogeneity BONE MARROW

## 05.12.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(Matrix)
library(genefilter)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(STutility)
library(ComplexHeatmap)
color <- rev(brewer.pal(11,"Spectral"))
source("./regress_out_function.R")

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

## Regress out, al final con mi funcion
st <- SCTransform(se.merged, assay = "RNA", new.assay.name = "SCT_novars", variable.features.n = 18000)
st <- SCTransform(se.merged, assay = "RNA", new.assay.name = "SCT_novars_sample", vars.to.regress = "name")
#ESTE
st <- regress_out(se.merged, "name")

## Aggregate counts by community en este caso con el regress por name
Seurat::Idents(object = st) <- st@meta.data[["community"]]
cts <- AggregateExpression(st, 
                           group.by = c("community"),
                           assays = 'regress',
                           #assays = 'SCT_novars_sample',
                           slot = "count",
                           return.seurat = FALSE)

aggregated_counts <- cts$regress
head(aggregated_counts)
row_totals <- rowSums(aggregated_counts)
filtered_data <- aggregated_counts[row_totals >= 10, ]
row_totals <- rowSums(filtered_data)
quantile_99 <- quantile(filtered_data, probs = 0.99)
#filtered_data <- filtered_data[row_totals <= 500, ]
head(filtered_data)

##for regress out
final_mat <-aggregated_counts

## Normalize by library
library_sizes <- colSums(final_mat)
scaling_factor <- mean(library_sizes)
normalized_matrix <- apply(final_mat, 2, library_sizes, FUN="/") * scaling_factor

final_mat <- normalized_matrix
head(final_mat)

## top var genes
### Estandar error
se <- function(x) sd(x)/sqrt(length(x))
t_final_mat <- Matrix::t(final_mat)
se_mat <- apply(t_final_mat, 2, se)
se_mat <- se_mat[order(se_mat, decreasing = T)]
se_mat_top100 <- names(sort(se_mat, decreasing = TRUE))[1:100]

## PCA
pca<- prcomp(t(final_mat),center = TRUE, scale. = TRUE) 
# check the order of the samples are the same.
all.equal(rownames(pca$x), colnames(final_mat))
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

top <- names(sort(se_mat, decreasing = TRUE))[1:100]

## Ordenar por los genes
final_mat_se <- final_mat[top, ] #pseudobulk ordered most var genes
dim(final_mat_se)

final_mat_se <- NormalizeData(final_mat_se)
final_mat_se <- ScaleData(final_mat_se)


# Draw the heatmap
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

suppressMessages(
  Heatmap(final_mat_se, show_column_names = TRUE, 
          row_names_gp = gpar(fontsize=6), cluster_columns = FALSE,
          top_annotation = top_annotation)
)

## Chequear genes
#noRegressGenes <- top
nameRegressGenes <- se_mat_top100
#pcRegreesGenes <- se_mat_top100

novsname <- intersect(nameRegressGenes, noRegressGenes)
top5 <- top[1:50]
top5


## Plots
x <- "Ly6d"
m1 <- FeatureOverlay(M1, features = c(x),pt.size = 1.3, col=color)
m2 <- FeatureOverlay(M2, features = c(x),pt.size = 1.3, col=color)
m3 <- FeatureOverlay(M8, features = c(x),pt.size = 1.3, col=color)
m4 <- FeatureOverlay(M9, features = c(x),pt.size = 1.3, col=color)

grid.arrange(m1, m2, m3, m4, ncol = 2)


se.merged <- NormalizeData(st)
se.merged <- ScaleData(st)
VlnPlot(se.merged, features = x, group.by = "community", assay = "SCT")



