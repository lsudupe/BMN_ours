## SCRIPT: Cluster 4vs5vs6 pseudobulk all samples Bone Marrow project

## 01.05.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(edgeR)
library(dplyr)
library(data.table)
library(sva)
library(STutility)
library(SingleCellExperiment)
library(Matrix.utils)

#Data---------------------------------
all <- readRDS("./objects/heterogeneity/se_hierarchical.rds")
b <- SetIdent(all, value = all@meta.data[["clustering"]])
b <- SubsetSTData(b, idents = c("5", "6","7"))

b <- SetIdent(b, value = b@meta.data[["name"]])
subset <- SubsetSTData(b, idents = c("M1_fem_1C", "M2_F_2B", "M8_F2_1C", "M9_F2_1C"))

b <- SetIdent(subset, value = subset@meta.data[["clustering"]])

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- b@assays$RNA@counts 
numbermorezero <- apply(counts,1,function(x){sum(x>1)})
counts <- counts[numbermorezero>2,]  

metadata <- b@meta.data

metadata$name[metadata$name == "M1_fem_1C"] <- "M1"
metadata$name[metadata$name == "M2_F_2B"] <- "M2"
metadata$name[metadata$name == "M8_F2_1C"] <- "M8"
metadata$name[metadata$name == "M9_F2_1C"] <- "M9"


# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(b@active.ident)
metadata$sample_id <- factor(metadata$name)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)
## Check the assays present
assays(sce)
## Check the counts matrix
dim(counts(sce))
counts(sce)[1:6, 1:6]
dim(colData(sce))
head(colData(sce))

###Preparing the single-cell dataset for pseudobulk analysis
# Extract unique names of clusters (= levels of cluster_id factor variable)
cluster_names <- levels(colData(sce)$cluster_id)
cluster_names

# Total number of clusters
length(cluster_names)

# Extract unique names of samples (= levels of sample_id factor variable)
sample_names <- levels(colData(sce)$sample_id)
sample_names

# Total number of samples
length(sample_names)

##Aggregating counts to the sample level for each cluster
# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("cluster_id", "sample_id")]
head(groups)

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum")

# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]

############## Splitting the counts matrix by cell type ##############
# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
aggr_counts[1:6, 1:6]

# check
df <- as.data.frame(aggr_counts)

## Exploring structure of function output (list
tstrsplit(colnames(aggr_counts), "_") %>% str()

## Comparing the first 10 elements of our input and output strings
head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)


#######create a list for each cluster
# As a reminder, we stored our cell types in a vector called cluster_names
cluster_names

# Loop over all cell types to extract corresponding counts, and store information in a list
## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

# Explore the different components of the list
str(counts_ls)
#######create a list for each cluster FIN

#######Generating matching metadata at the sample-level
# Reminder: explore structure of metadata
head(colData(sce))
# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(cluster_id, sample_id)

dim(metadata)
head(metadata)

# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]
dim(metadata)
head(metadata)
# Rename rows
rownames(metadata) <- apply(metadata, 1, function(row) { paste(row[1], row[2], sep = "_") })
head(metadata)

# Number of spots per sample and cluster
t <- table(colData(sce)$sample_id,
           colData(sce)$cluster_id)
#t[1:6, 1:6]

###we will append this cell count information to our generic metadata table
# Creating metadata list

## Initiate empty list
metadata_ls <- list()
for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}

# Explore the different components of the list
str(metadata_ls)

####Differential gene expression with DESeq2#####
##Creating a DESeq2 object
# Select cell type of interest
cluster_names

# Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))

# cluster 5
idx5 <- which(names(counts_ls) == "5")
cluster_counts5 <- counts_ls[[idx5]]
cluster_metadata5 <- metadata_ls[[idx5]]
# cluster 6
idx6 <- which(names(counts_ls) == "6")
cluster_counts6 <- counts_ls[[idx6]]
cluster_metadata6 <- metadata_ls[[idx6]]
# cluster 7
idx7 <- which(names(counts_ls) == "7")
cluster_counts7 <- counts_ls[[idx7]]
cluster_metadata7 <- metadata_ls[[idx7]]

# Check contents of extracted objects for cluster 5
head(cluster_metadata5)
# Check contents of extracted objects for cluster 6
head(cluster_metadata6)
# Check contents of extracted objects for cluster 7
head(cluster_metadata7)

# concatenate counts data for all three clusters
cluster_counts_all <- do.call(cbind, list(cluster_counts5, cluster_counts6,cluster_counts7))
# concatenate metadata data for all three clusters
cluster_metadata_all <- rbind(cluster_metadata5, cluster_metadata6,cluster_metadata7)

# Check matching of matrix columns and metadata rows
all(colnames(cluster_counts_all) == rownames(cluster_metadata_all))

# Create DESeq2 object        
dds <- DESeqDataSetFromMatrix(cluster_counts_all, 
                              colData = cluster_metadata_all, 
                              design = ~ cluster_id)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, ntop = 500, intgroup = "cluster_id")

#Hierarchical clustering
# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata_all[, c("cluster_id"), drop=F])

#Running DESeq2
# Run DESeq2 differential expression analysis
dds <- DESeq(dds)
# Plot dispersion estimates
plotDispEsts(dds)

######Exploring DE results
# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, 
               name = "cluster_id_6_vs_5",
               alpha = 0.05)

res <- results(dds)



#####Table of significant for all genes
# Turn the DESeq2 results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)


# Check results output
res_tbl 
# Set thresholds
padj_cutoff <- 0.05
# Subset the significant results
sig_res <- dplyr::filter(res_tbl, pvalue < padj_cutoff) %>%
  dplyr::arrange(padj)
# Check significant genes output
sig_res

#####Heatmap
# Heatmap
## Extract normalized counts from dds object
normalized_counts <- counts(dds, normalized = TRUE)
## Extract normalized counts for significant genes only
sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]

## Set a color-blind friendly palette
heat_colors <- rev(brewer.pal(11, "PuOr"))

my_col_order <- c("4_M1", "4_M2", "4_M8" ,"4_M9" ,"5_M1" ,"5_M2" ,"5_M8", "5_M9", "6_M2", "6_M8")

## Run pheatmap using the metadata data frame for the annotation
pdf("./results/DE/st/pseudo/pseudo_normal.pdf")
pheatmap(sig_counts, 
         color = heat_colors, 
         cluster_rows = TRUE, 
         show_rownames = FALSE,
         annotation = cluster_metadata_all[, c("sample_id", "cluster_id")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)  
dev.off()

#####normal DE with covariate

# Differential expression analysis between group1 and group2
diff_exp_results_12 <- FindMarkers(b, ident.1 = '5', ident.2 = '6', latent.vars = "MM_MIC",
                                   min.cells.feature = 3,
                                   min.cells.group = 3)

# Differential expression analysis between group1 and group3
diff_exp_results_13 <- FindMarkers(b, ident.1 = '5', ident.2 = '7', latent.vars = "MM_MIC",
                                   min.cells.feature = 3,
                                   min.cells.group = 3)

# Differential expression analysis between group2 and group3
diff_exp_results_23 <- FindMarkers(b, ident.1 = '6', ident.2 = '7', latent.vars = "MM_MIC",
                                   min.cells.feature = 3,
                                   min.cells.group = 3)

# For each result object, filter based on adjusted p-value and log fold change
# For each result object, filter based on adjusted p-value and log fold change
significant_results_12 <- diff_exp_results_12[diff_exp_results_12$p_val_adj < 0.05 & abs(diff_exp_results_12$avg_log2FC) > 1,]
significant_results_13 <- diff_exp_results_13[diff_exp_results_13$p_val_adj < 0.05 & abs(diff_exp_results_13$avg_log2FC) > 1,]
significant_results_23 <- diff_exp_results_23[diff_exp_results_23$p_val_adj < 0.05 & abs(diff_exp_results_23$avg_log2FC) > 1,]


# Add a new column to the results indicating whether each gene is highly significant
diff_exp_results_12$highly_sig <- diff_exp_results_12$p_val_adj < 0.05 & abs(diff_exp_results_12$avg_log2FC) > 1

library(ggrepel)
library(gggenes)
# Volcano plot with different colors for highly significant genes
pdf("./results/DE/st/volcano5vs6.pdf")
ggplot(diff_exp_results_12, aes(x = avg_log2FC, y = -log10(p_val_adj), color = highly_sig)) +
  geom_point() +
  geom_text_repel(data = subset(diff_exp_results_12, highly_sig), 
                  aes(label = rownames(subset(diff_exp_results_12, highly_sig))), 
                  size = 3) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic() +
  xlab("Log2 fold change") +
  ylab("-Log10 adjusted p-value") +
  ggtitle("Volcano plot: Group5 vs Group6")
dev.off()

diff_exp_results_13$highly_sig <- diff_exp_results_13$p_val_adj < 0.05 & abs(diff_exp_results_13$avg_log2FC) > 1

pdf("./results/DE/st/volcano5vs7.pdf")
ggplot(diff_exp_results_13, aes(x = avg_log2FC, y = -log10(p_val_adj), color = highly_sig)) +
  geom_point() +
  geom_text_repel(data = subset(diff_exp_results_13, highly_sig), 
                  aes(label = rownames(subset(diff_exp_results_13, highly_sig))), 
                  size = 3) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic() +
  xlab("Log2 fold change") +
  ylab("-Log10 adjusted p-value") +
  ggtitle("Volcano plot: Group5 vs Group7")
dev.off()

diff_exp_results_23$highly_sig <- diff_exp_results_23$p_val_adj < 0.05 & abs(diff_exp_results_23$avg_log2FC) > 1

pdf("./results/DE/st/volcano6vs7.pdf")
ggplot(diff_exp_results_23, aes(x = avg_log2FC, y = -log10(p_val_adj), color = highly_sig)) +
  geom_point() +
  geom_text_repel(data = subset(diff_exp_results_23, highly_sig), 
                  aes(label = rownames(subset(diff_exp_results_23, highly_sig))), 
                  size = 3) +
  scale_color_manual(values = c("black", "red")) +
  theme_classic() +
  xlab("Log2 fold change") +
  ylab("-Log10 adjusted p-value") +
  ggtitle("Volcano plot: Group6 vs Group7")
dev.off()

# Extract the names of the "highly significant" genes
highly_sig_genes_12 <- rownames(subset(diff_exp_results_12, highly_sig))
highly_sig_genes_13 <- rownames(subset(diff_exp_results_13, highly_sig))
highly_sig_genes_23 <- rownames(subset(diff_exp_results_23, highly_sig))
# Print the names
print(highly_sig_genes_12)
print(highly_sig_genes_13)
print(highly_sig_genes_23)

highly_sig_genes <- c(highly_sig_genes_12, highly_sig_genes_23)

##save data
# Write combined_lines to a .txt file
writeLines(highly_sig_genes_12, "./results/DE/st/5vs6.txt")
writeLines(highly_sig_genes_13, "./results/DE/st/5vs7.txt")
writeLines(highly_sig_genes_23, "./results/DE/st/6vs7.txt")

#spatial plots
color <- brewer.pal(11,"Spectral")
color <- rev(color)

pdf(file.path("./results/DE/st/Ighj4_all.pdf"))
FeatureOverlay(all, features = "Ighj4",sampleids = 1:6, ncols = 2,pt.size = 0.7,cols = color)
dev.off()

gene_expression <- b@assays$RNA@counts[highly_sig_genes, ]
# Scale the data
scaled_data <- scale(t(gene_expression))

# Create the heatmap
pdf("./results/DE/st/significant_heatmap.pdf")
DoHeatmap(b, features = highly_sig_genes) +
  theme_classic() +
  xlab("Groups") +
  ylab("Top Differentially Expressed Genes") +
  ggtitle("Heatmap")
dev.off()

#####normal DE with covariate FIN

####GO analysis
library(clusterProfiler)
library(org.Mm.eg.db)

##foldchange decreciente!!!
##no filtrar por pvalue
##lista de patways enriqucidos

genes_12 <- diff_exp_results_12[diff_exp_results_12$p_val_adj < 0.05 & abs(diff_exp_results_12$avg_log2FC) > 1.5,]
genes_13 <- diff_exp_results_13[diff_exp_results_13$p_val_adj < 0.05 & abs(diff_exp_results_13$avg_log2FC) > 1.5,]
genes_23 <- diff_exp_results_23[diff_exp_results_23$p_val_adj < 0.05 & abs(diff_exp_results_23$avg_log2FC) > 1.5,]


genes_12 = bitr(rownames(genes_12), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
genes_13 = bitr(rownames(genes_13), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
genes_23 = bitr(rownames(genes_23), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


genelist <- list(genes_12,
                 genes_13,
                 genes_23)
                 
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
ck <- setReadable(GOclusterplot, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
                 
pdf("./results/DE/st/se_all_clusters_dotplot_rna.pdf")
print(dotplot(GOclusterplot))
dev.off()
                 
####GO analysis FIN

library(clusterProfiler)
library(org.Mm.eg.db)

group5vs6 <- diff_exp_results_12
group5vs6 <- group5vs6[order(group5vs6$avg_log2FC, decreasing = TRUE),]
group5vs6_ = bitr(rownames(group5vs6), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
group5vs7 <-diff_exp_results_13
group5vs7 <- group5vs7[order(group5vs7$avg_log2FC, decreasing = TRUE),]
group5vs7_ = bitr(rownames(group5vs7), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
group6vs7 <-diff_exp_results_23
group6vs7 <- group6vs7[order(group6vs7$avg_log2FC, decreasing = TRUE),]
group6vs7_ = bitr(rownames(group6vs7), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


# list of Seurat DE results
listSeuratDE <- list("group5vs6" = group5vs6_$ENTREZID,
                 "group5vs7" = group5vs7_$ENTREZID,
                 "group6vs7" = group6vs7_$ENTREZID)

# Assuming 'df' is one of your differential expression result data frames
group6vs7_prueba <- df[order(df$log2FoldChange, decreasing = TRUE),]

# Create a named numeric vector
geneList <- df$log2FoldChange
names(geneList) <- df$ENTREZID


GOclusterplot <- compareCluster(geneCluster = listSeuratDE, fun = "enrichGO", OrgDb = "org.Mm.eg.db",)
GOenrichment <- compareCluster(listSeuratDE, fun = "gseGO",OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
KEGGenrichment <- compareCluster(listSeuratDE, fun = "gseKEGG", organism = "mmu")


listSeuratDE <- list(group5vs6, group5vs7, group6vs7)  
listDiffExprFilter <- list()
orderGenes <- list()

for (i in 1:length(listSeuratDE)){
  # Add ENTREZ id
  listSeuratDE[[i]]$ENTREZ <- mapIds(x=org.Mm.eg.db, column = "ENTREZID", 
                                     key = rownames(listSeuratDE[[i]]), keytype = "SYMBOL")
  
  # Filter and order by avg_log2FC
  listDiffExprFilter[[i]] <- listSeuratDE[[i]][order(listSeuratDE[[i]]$avg_log2FC, decreasing = T),]
  
  # Filter out rows with NA in avg_log2FC or ENTREZ
  listDiffExprFilter[[i]] <- listDiffExprFilter[[i]][!(is.na(listDiffExprFilter[[i]]$avg_log2FC)) | !is.na(listDiffExprFilter[[i]]$ENTREZ),] 
  
  # Store the ordered avg_log2FC values with their corresponding ENTREZ ids
  orderGenes[[i]] <- listDiffExprFilter[[i]]$avg_log2FC
  names(orderGenes[[i]]) <- listDiffExprFilter[[i]]$ENTREZ
  
}

names(orderGenes) <- paste0("Cluster", 1:length(orderGenes))

GOenrichment <- compareCluster(orderGenes, fun = "gseGO",OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
ck <- setReadable(GOenrichment, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
dotplot(GOenrichment, showCategory=5) + ggtitle("dotplot for GSEA")
dotplot(ck)

KEGGenrichment <- compareCluster(orderGenes, fun = "gseKEGG", organism = "mmu")

source("./PlotTools_Dani.R")
## plot_Dani is a dot plot showing only the top pathways for each comparison
plot_Dani(GOenrichment)
plot_Dani(KEGGenrichment)
### plot_heatmpaGSEA_Dani is a heatmap plot showing all the enriched pathway for each comparison
plot_heatmapGSEA_Dani(GOenrichment)
plot_heatmapGSEA_Dani(KEGGenrichment)

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

head(listSeuratDE[[1]])

