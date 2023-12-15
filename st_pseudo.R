## SCRIPT: Heterogeneity BONE MARROW

## 05.12.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(Matrix)
library(sva)
library(genefilter)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(STutility)
library(ComplexHeatmap)
library(DESeq2)
color <- rev(brewer.pal(11,"Spectral"))
#source("./regress_out_function.R")

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

## quitar inmunos
#immunoglobulin_genes <- grep("^Igh", rownames(se.merged), value = TRUE)
# Remove these genes from the matrix
#se.merged <- se.merged[!rownames(se.merged) %in% immunoglobulin_genes, ]

## Aggregate counts by community en este caso con el regress por name
Seurat::Idents(object = se.merged) <- se.merged@meta.data[["community"]]
cts <- AggregateExpression(se.merged, 
                           group.by = c("community"),
                           assays = 'RNA',
                           slot = "count",
                           return.seurat = FALSE)

aggregated_counts <- cts$RNA
head(aggregated_counts)

## regress out sample with Combat
batch<- c("S1", "S1", "S2", "S2","S2","S2","S2" ,"S8","S8","S8","S9","S9")
modules <- c("M1", "M1","M2","M2","M2","M2","M2","M8","M8","M9","M9","M9")
corrected_matrix <- ComBat_seq(counts=aggregated_counts, batch = batch, group = NULL)

head(corrected_matrix)
corrected_matrix <- aggregated_counts

colData_df <- DataFrame(modules)
dds_matrix <- DESeqDataSetFromMatrix(countData = corrected_matrix, colData = colData_df, design = ~modules)
dds_matrix <- DESeq(dds_matrix)

## Extract normalize counts
final_mat <- counts(dds_matrix, normalized = T)
head(final_mat)

## Estandar error
se <- function(x) sd(x)/sqrt(length(x))
t_final_mat <- Matrix::t(final_mat)
se_mat <- apply(t_final_mat, 2, se)
se_mat <- se_mat[order(se_mat, decreasing = T)]
hist(se_mat, main="Histogram of se_mat", xlab="standard error values")
log_se_mat <- log1p(se_mat) # log1p is used to avoid log(0) which is undefined
hist(log_se_mat, main="Log-transformed Histogram of se_mat", xlab="Log-transformed  standard error values")

## select top
# Calculate the 99th quantile on the log-transformed data
quantile_99 <- quantile(log_se_mat, probs = 0.99)
threshold_value <- quantile_99[1]
# To get the threshold value on the original scale of se_mat
original_threshold_value <- exp(threshold_value) - 1
selected_genes <- se_mat[se_mat > original_threshold_value]

## Ordenar por los genes
top <- names(sort(selected_genes, decreasing = TRUE))#[1:100]
final_mat_se <- final_mat[top, ] #pseudobulk ordered most var genes
dim(final_mat_se)

## ana heatmap 
matrix <- final_mat[top, ]
ordered_cols <- c("S1_M1", "S1_M2", "S2_M1", "S2_M2", "S2_M3", "S2_M4", "S2_M5", 
                  "S8_M1", "S8_M2", "S9_M1", "S9_M2", "S9_M3")
matrix_ordered <- matrix[, ordered_cols]

pheatmap(matrix_ordered, scale = "row", fontsize = 4)

## pca
pca<- prcomp(t(matrix),center = FALSE, scale. = FALSE) 
PC1_and_PC2<- data.frame(PC1=pca$x[,1], PC2= pca$x[,2])
PC1_and_PC2$sample <- c("M1","M1", "M2", "M2","M2","M2","M2", "M8","M8","M9","M9","M9")
head(PC1_and_PC2)

ggplot(PC1_and_PC2, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color = sample)) +
  theme_bw(base_size = 14) 

ggplot(PC1_and_PC2, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color = rownames(PC1_and_PC2))) +
  theme_bw(base_size = 14) 
 
## plots spatial y violin
top5 <- top[1:100]
top5

top_batcheffect_name <- top
top_NObatcheffect_name <- top

genesComun <- intersect(top_batcheffect_name, top_NObatcheffect_name)

x <- "Mzb1"
m1 <- FeatureOverlay(M1, features = c(x),pt.size = 1.3, col=color)
m2 <- FeatureOverlay(M2, features = c(x),pt.size = 1.3, col=color)
m3 <- FeatureOverlay(M8, features = c(x),pt.size = 1.3, col=color)
m4 <- FeatureOverlay(M9, features = c(x),pt.size = 1.3, col=color)

grid.arrange(m1, m2, m3, m4, ncol = 2)

"Xbp1" %in% top_NObatcheffect_name
"Mzb1" %in% top_NObatcheffect_name

VlnPlot(se.merged, features = x, group.by = "community", assay = "SCT")

FeatureOverlay(se.merged, features = c("community"), 
               sampleids = 1:4, ncols = 2, 
               pt.size = 1.3)

FeatureOverlay(se.merged, features = "community")


p1 <- ST.FeaturePlot(se.merged, features = "community", indices = 1,  pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
p2 <- ST.FeaturePlot(se.merged, features = "community", indices = 2, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
p3 <- ST.FeaturePlot(se.merged, features = "community", indices = 3, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
p4 <- ST.FeaturePlot(se.merged, features = "community", indices = 4, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
cowplot::plot_grid(p1, p2, p3, p4, ncol = 4)

ST.FeaturePlot(M1, features = "community", split.labels = T, pt.size = 1.5) & theme(plot.title = element_blank(), strip.text = element_blank())
ST.FeaturePlot(M2, features = "community", split.labels = T, pt.size = 1.5) & theme(plot.title = element_blank(), strip.text = element_blank())


FeatureOverlay(M9, features = "community")






