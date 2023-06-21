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

#############REGRESSOUT
met <- b@meta.data
df_grouped <- met %>% 
  group_by(name, clustering) %>% 
  summarize(mean_MM_MIC = mean(MM_MIC)) %>% 
  ungroup()
# pivot the dataframe to get name as row names, clustering as columns, and MM_MIC as values
df_pivot <- pivot_wider(df_grouped, 
                        names_from = clustering, 
                        values_from = mean_MM_MIC, 
                        names_prefix = "cluster_",
                        values_fill = 0) %>% 
  column_to_rownames("name")
# sort the row names alphabetically
df_pivot <- df_pivot[order(rownames(df_pivot)),]
df_pivot <- data.frame(lapply(df_pivot, function(x) x*100))
# display the resulting dataframe
print(df_pivot)


##cluster5##
# extract the values of column of interest from df_pivot and add them as a last row in df
df_5 <- as.data.frame(cluster_counts5)
matrix_t <- t(df_5)
matrix_t <- as.data.frame(matrix_t)
matrix_t["plasma_value"] <- df_pivot$cluster_5
matrix_final <- as.data.frame(t(matrix_t))

###REGRESS OUT
data <- matrix_final
# Get the number of genes and spots
genes <- nrow(data) - 1 # Subtract 1 to exclude the last row with plasma cell percentage values
spots <- ncol(data)
# Convert the data frame to a matrix
data_matrix <- as.matrix(data)
# Initialize the residuals matrix with the same dimensions as the data matrix (without the last row)
residuals_matrix <- matrix(0, nrow = genes, ncol = spots)
# Loop through each gene
for (i in 1:genes) {
  # Create a linear model to regress out the plasma cell percentage
  model <- lm(data_matrix[i, 1:spots] ~ data_matrix[genes + 1, 1:spots])
  # Save the residuals in the corresponding row of the residuals_matrix
  residuals_matrix[i, ] <- model$residuals
}
# Replace the original data_matrix with the residuals
data_matrix[1:genes, 1:spots] <- residuals_matrix
# Convert the matrix back to a data frame
data <- as.data.frame(data_matrix)
data <- data[1:(nrow(data)-1),]

#for ana
cluster_counts5 <- as.matrix(matrix_final)
#for ana fin 

cluster_counts5 <- as.matrix(data)
cluster_counts5 <- Matrix(cluster_counts5, sparse = FALSE)

##cluster6##
# extract the values of column of interest from df_pivot and add them as a last row in df
df_6 <- as.data.frame(cluster_counts6)
matrix_t <- t(df_6)
matrix_t <- as.data.frame(matrix_t)
matrix_t["plasma_value"] <- df_pivot$cluster_6[2:3]
matrix_final <- as.data.frame(t(matrix_t))

###REGRESS OUT
data <- matrix_final
# Get the number of genes and spots
genes <- nrow(data) - 1 # Subtract 1 to exclude the last row with plasma cell percentage values
spots <- ncol(data)
# Convert the data frame to a matrix
data_matrix <- as.matrix(data)
# Initialize the residuals matrix with the same dimensions as the data matrix (without the last row)
residuals_matrix <- matrix(0, nrow = genes, ncol = spots)
# Loop through each gene
for (i in 1:genes) {
  # Create a linear model to regress out the plasma cell percentage
  model <- lm(data_matrix[i, 1:spots] ~ data_matrix[genes + 1, 1:spots])
  # Save the residuals in the corresponding row of the residuals_matrix
  residuals_matrix[i, ] <- model$residuals
}
# Replace the original data_matrix with the residuals
data_matrix[1:genes, 1:spots] <- residuals_matrix
# Convert the matrix back to a data frame
data <- as.data.frame(data_matrix)
data <- data[1:(nrow(data)-1),]

#for ana
cluster_counts6 <- as.matrix(matrix_final)
#for ana fin 

cluster_counts6 <- as.matrix(data)
cluster_counts6 <- Matrix(cluster_counts6, sparse = FALSE)

##cluster7##
# extract the values of column of interest from df_pivot and add them as a last row in df
df_7 <- as.data.frame(cluster_counts7)
matrix_t <- t(df_7)
matrix_t <- as.data.frame(matrix_t)
matrix_t["plasma_value"] <- df_pivot$cluster_7[1:3]
matrix_final <- as.data.frame(t(matrix_t))

###REGRESS OUT
data <- matrix_final
# Get the number of genes and spots
genes <- nrow(data) - 1 # Subtract 1 to exclude the last row with plasma cell percentage values
spots <- ncol(data)
# Convert the data frame to a matrix
data_matrix <- as.matrix(data)
# Initialize the residuals matrix with the same dimensions as the data matrix (without the last row)
residuals_matrix <- matrix(0, nrow = genes, ncol = spots)
# Loop through each gene
for (i in 1:genes) {
  # Create a linear model to regress out the plasma cell percentage
  model <- lm(data_matrix[i, 1:spots] ~ data_matrix[genes + 1, 1:spots])
  # Save the residuals in the corresponding row of the residuals_matrix
  residuals_matrix[i, ] <- model$residuals
}
# Replace the original data_matrix with the residuals
data_matrix[1:genes, 1:spots] <- residuals_matrix
# Convert the matrix back to a data frame
data <- as.data.frame(data_matrix)
data <- data[1:(nrow(data)-1),]

#for ana
cluster_counts7 <- as.matrix(matrix_final)
#for ana fin 

cluster_counts7 <- as.matrix(data)
cluster_counts7 <- Matrix(cluster_counts7, sparse = FALSE)

#############REGRESSOUT FIN
# concatenate counts data for all three clusters
cluster_counts_all <- do.call(cbind, list(cluster_counts5, cluster_counts6, cluster_counts7))

#save data for ana
write.csv(cluster_counts_all, "./counts.csv", row.names=TRUE)
write.csv(cluster_metadata_all, "./meta.csv", row.names=TRUE)

# scale the values in the matrix to be from 0 to 1
scaled_mat <- apply(cluster_counts_all, 2, function(x) (x - min(x)) / (max(x) - min(x)))

mat_integer <- as.matrix(as.integer(scaled_mat))

# Create DESeq2 object        
dds <- DESeqDataSetFromMatrix(scaled_mat, 
                              colData = cluster_metadata_all, 
                              design = ~ cluster_id)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, ntop = 500, intgroup = "cluster_id")
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
               name = "cluster_id_5_vs_4",
               alpha = 0.05)

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
res <- lfcShrink(dds, 
                 coef = "group_id_stim_vs_ctrl",
                 res=res,
                 type = "apeglm")

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
pdf("./results/DE/st/pseudo/pseudo_regress.pdf")
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





