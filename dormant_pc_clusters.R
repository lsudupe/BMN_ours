## SCRIPT: Dormant genes individual plot  BM project

## 01.06.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(dplyr)
library(UCell)
library(tidyr)

#Data---------------------------------
objects <- readRDS("./objects/sp/regress_out/list_regressout_ST.rds")
all <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

##individual
M1_s <- readRDS("./objects/pc_clusters/M1_s_dormant.rds")
M2_s <- readRDS("./objects/pc_clusters/M2_s_dormant.rds")
M8_s <- readRDS("./objects/pc_clusters/M8_s_dormant.rds")
M9_s <- readRDS("./objects/pc_clusters/M9_s_dormant.rds")


dormant <- c("C1qa", "Aif1" ,"Axl" ,"II18bp", "Glul", "Mpeg1", "H2-Eb1",
             "Nr1h3" ,"Lilrb4" ,"Sdc3" ,"Fcgr3" ,"Oas1a", "Ifit2", "Ly6a",
             "Ly86" ,"Bcl2a1b", "Cd19" ,"Gm4955", "Hebp1", "Sirpa" ,"Slfn2" ,"Tnfaip2" ,
             "Vpreb3", "ligp1", "Oasl2","am213b" ,"Slc43a2", "Ifi44", "Gng10" ,"Cxcl16" ,"Ptprcap" ,
             "Anxa2" ,"Rgs2" ,"Tmed3" ,"Igll1", "Hpgd", "Glipr1", "Cd4" ,"Cd84" ,"Gbp2", "AB124611", "Slc44a2" ,
             "Samd9l" ,"Oas1g" ,"Fcgr1", "Pla2g15", "Tifa" ,"Pmp22", "Abcc3" ,"S100a10")

jose_genes <- c ("Cd44", "Cd81", "Flna", "Mki67", "Pcna", "Xbp1")

#M1
Idents(object = M1_s) <- M1_s@meta.data[["pc_clusters"]]
M1_subset <- SubsetSTData(object = M1_s, ident = c("M1_group_1_cluster5","M1_group_1_cluster6", "M1_group_1_cluster7"))

#M2
Idents(object = M2_s) <- M2_s@meta.data[["pc_clusters"]]
M2_subset <- SubsetSTData(object = M2_s, ident = c("M2_group_1_cluster5","M2_group_1_cluster6", "M2_group_1_cluster7",
                                                   "M2_group_2_cluster5", "M2_group_2_cluster6", "M2_group_2_cluster7",
                                                   "M2_group_3_cluster5","M2_group_3_cluster6", "M2_group_3_cluster7",
                                                   "M2_group_4_cluster5","M2_group_4_cluster6", "M2_group_4_cluster7"
))

#M8
Idents(object = M8_s) <- M8_s@meta.data[["pc_clusters"]]
M8_subset <- SubsetSTData(object = M8_s, ident = c("M8_group_1_cluster5","M8_group_1_cluster6", "M8_group_1_cluster7",
                                                   "M8_group_2_cluster5", "M8_group_2_cluster6", "M8_group_2_cluster7"))


#M9
Idents(object = M9_s) <- M9_s@meta.data[["pc_clusters"]]
M9_subset <- SubsetSTData(object = M9_s, ident = c("M9_group_1_cluster6",
                                                   "M9_group_2_cluster5", "M9_group_2_cluster6",
                                                   "M9_group_3_cluster6"))

#######TWO LOOPS FOR TWO LISTS
normal <- c(M1_s,M2_s,M8_s,M9_s)
names(normal) <- c("M1_s","M2_s","M8_s","M9_s")

for (i in 1:length(normal)){
  a <- normal[[i]]
  b <- names(normal[i])
  
  color <- rainbow(30)
  ######
  # Get the expression data
  data <- FetchData(a, vars = jose_genes)
  
  # Add the cluster data to the data frame
  data$cluster <- a@meta.data$labels
  
  # Transform the data into a long format
  df <- gather(data, gene, expression, -cluster)
  
  # Remove genes with zero expression
  df <- df %>% filter(expression > 0)
  
  # Plot
  plot <- ggplot(df, aes(x = gene, y = expression, fill = cluster)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Gene", y = "Gene Expression") +
    coord_flip()
  
  pdf(paste("./results/pc_clusters/jose_genes/all/", b,"_j_a_genes_labels.pdf",sep=""), width = 10, height = 7)
  print(plot)
  dev.off()
  
  ######
  #p <- FeatureOverlay(a, features = c("pc_clusters"), ncols = 1, pt.size = 1.5,cols = color)
  
  #pdf(paste("./results/pc_clusters/all/", b,"_clusters567.pdf",sep=""))
  #print(p)
  #dev.off()
  
  # Boxplot
  df <- a@meta.data
  df_selected <- df %>% select(name, labels, signature_1_dormant, labels)
  
  plot <- ggplot(df_selected, aes(x = labels, y = signature_1_dormant, fill = labels)) +
    geom_boxplot() +
    labs(title = "Dormant distributions with respect to labels",
         x = "Groups", 
         y = "Signature Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/jose_genes/all/", b,"_dormant_labels.pdf",sep=""), width = 10, height = 7)
  print(plot)
  dev.off()
  
  ##dormant signature in clusters boxplot 
  
  df_selected <- df %>% select(name, clustering, signature_1_dormant, clustering)
  
  plot <- ggplot(df_selected, aes(x = clustering, y = signature_1_dormant, fill = clustering)) +
    geom_boxplot() +
    labs(title = "Dormant distributions with respect to clusters",
         x = "Groups", 
         y = "Signature Dormant", 
         fill = "clusters") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/clusters/", b,"_dormant_clusters.pdf",sep=""), width = 10, height = 7)
  print(plot)
  dev.off()
  
  #residuals_dormant_ucellscore
  
  df_selected <- df %>% select(name, clustering, residuals_dormant_ucellscore, clustering)
  
  plot <- ggplot(df_selected, aes(x = clustering, y = residuals_dormant_ucellscore, fill = clustering)) +
    geom_boxplot() +
    labs(title = "Dormant residuals distributions with respect to clusters",
         x = "Groups", 
         y = "Signature Dormant", 
         fill = "clusters") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/clusters/", b,"_dormant_residuals_clusters.pdf",sep=""), width = 10, height = 7)
  print(plot)
  dev.off()
  
}



subsets <- c(M1_subset,M2_subset, M8_subset, M9_subset)
names(subsets) <- c("M1_subset","M2_subset", "M8_subset", "M9_subset")

for (i in 1:length(subsets)){
  a <- subsets[[i]]
  b <- names(subsets[i])
  
  color <- rainbow(30)
  
  ######
  # Get the expression data
  data <- FetchData(a, vars = jose_genes)
  
  # Add the cluster data to the data frame
  data$cluster <- a@meta.data$labels
  
  # Transform the data into a long format
  df <- gather(data, gene, expression, -cluster)
  
  # Remove genes with zero expression
  df <- df %>% filter(expression > 0)
  
  # Plot
  plot <- ggplot(df, aes(x = gene, y = expression, fill = cluster)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Gene", y = "Gene Expression") +
    coord_flip()
  
  pdf(paste("./results/pc_clusters/jose_genes/subset/", b,"_j_a_genes_labels.pdf",sep=""), width = 10, height = 7)
  print(plot)
  dev.off()
  
  ######
  #p <- FeatureOverlay(a, features = c("pc_clusters"), ncols = 1, pt.size = 1.5,cols = color)
  
  #pdf(paste("./results/pc_clusters/all/", b,"_clusters567.pdf",sep=""))
  #print(p)
  #dev.off()
  
  df <- a@meta.data
  df_selected <- df %>% select(name, labels, signature_1_dormant, labels)
  
  plot <- ggplot(df_selected, aes(x = labels, y = signature_1_dormant, fill = labels)) +
    geom_boxplot() +
    labs(title = "Dormant distributions with respect to labels",
         x = "Groups", 
         y = "Signature Dormant", 
         fill = "Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +  # Adjust text size and angle here
    facet_wrap(~name)  # Separate plots by 'name'
  
  pdf(paste("./results/pc_clusters/jose_genes/subset/", b,"_dormant_labels.pdf",sep=""), width = 10, height = 7)
  print(plot)
  dev.off()
  
}

###merge data
se_merged <- MergeSTData(M1_s, y = c(M2_s, M8_s, M9_s), 
                    add.spot.ids = c("M1_s", "M2_s", "M8_s", "M9_s"), project = "BM")


###correlation part
meta <- se_merged@meta.data

colnames <- colnames(meta)


##1. Calculating correlation values
# Select relevant columns for correlation analysis
cols_of_interest <- c("Tcell", "Bcell", "MM_MIC", "Erythroblasts", "Monocytes", "Neutrophils",
                      "DC", "signature_1_dormant", "residuals_dormant_ucellscore")

# Create a subset of the dataframe with the selected columns
subset_df <- meta[, cols_of_interest]

# Calculate correlation matrix
correlation_matrix <- cor(subset_df)

# Write the correlation matrix to an Excel file
write.csv(correlation_matrix, file = "correlation_matrix.csv", row.names = TRUE)


#2: Adding gene expression values and calculating correlations
# Extract normalized expression values for selected genes
jose_genes <- c("Cd44", "Cd81", "Flna", "Mki67", "Pcna", "Xbp1")
gene_expression <- FetchData(se_merged, vars = jose_genes)

# Merge gene expression values with the existing dataframe
merged_df <- cbind(meta, gene_expression)

# Calculate correlation between gene expression values and groups
gene_correlation <- cor(merged_df[, c("Tcell", "Bcell", "MM_MIC", "Erythroblasts", "Monocytes", 
                                      "Neutrophils", "DC", "signature_1_dormant",
                                      "residuals_dormant_ucellscore")],
                        merged_df[, jose_genes])

# Write the gene correlation matrix to an Excel file
write.csv(gene_correlation, file = "gene_correlation_matrix.csv", row.names = TRUE)


#3: Plotting correlation matrices
# Read the correlation matrix CSV files
# Read the correlation matrix CSV files
correlation_matrix <- read.csv("correlation_matrix.csv", row.names = 1)
gene_correlation <- read.csv("gene_correlation_matrix.csv", row.names = 1)

# Convert correlation matrix values to numeric
correlation_matrix[] <- lapply(correlation_matrix, as.numeric)
gene_correlation[] <- lapply(gene_correlation, as.numeric)

correlation_matrix <- as.matrix(correlation_matrix)
gene_correlation <- as.matrix(gene_correlation)

# Plot correlation matrix heatmap
heatmap_plot <- heatmap(correlation_matrix, col = colorRampPalette(c("white", "blue")))
heatmap_plot <- heatmap(correlation_matrix, col = colorRampPalette(c("white", "blue"))(100))

# Plot gene correlation matrix heatmap
gene_heatmap_plot <- heatmap(gene_correlation, col = colorRampPalette(c("white", "blue")))

# Display the heatmaps
library(gplots)
heatmap.2(correlation_matrix, col = colorRampPalette(c("white", "blue")))

plot(heatmap_plot)
plot(gene_heatmap_plot)














