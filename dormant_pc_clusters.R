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

jose_genes <- c ("Cd44", "Cd81", "Flna", "Mki67", "Pcna", "Xbp1", "Cdc37")

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

# Assume 'seurat_obj' is your Seurat object
# Normalize data
seurat_obj <- NormalizeData(se_merged)

# Fetch the normalized expression data for each gene
jose_genes <- c("Cd44", "Cd81", "Flna", "Mki67", "Pcna", "Xbp1", "Cdc37")
gene_data <- FetchData(seurat_obj, vars = jose_genes)

df <- seurat_obj@meta.data
# Combine the data
df <- cbind(df, gene_data)

###########CORRELATION PLOTS AND EXCELL#########
install.packages(c('tidyverse', 'readxl', 'writexl'))
library(openxlsx)
library(tidyverse)

# Function to create an Excel file with separate sheets for each group within a category
create_excel <- function(df, category) {
  # Get unique groups within the category
  groups <- unique(df[[category]])
  
  # Create a new Excel workbook
  wb <- createWorkbook()
  
  # For each group, compute the correlation matrix and add it to the workbook as a new sheet
  for (group in groups) {
    # Subset the data for this group
    df_group <- df %>% filter((!!sym(category)) == group)
    
    # Compute the correlation matrix after removing columns with zero standard deviation
    correlation_matrix <- df_group %>% 
      select(Tcell:DC, signature_1_dormant, residuals_dormant_ucellscore, all_of(jose_genes)) %>% 
      select_if(function(x) sd(x, na.rm=TRUE) != 0) %>% 
      cor()
    
    # Create a new sheet in the workbook for this correlation matrix
    addWorksheet(wb, paste("Group", group))
    writeData(wb, paste("Group", group), correlation_matrix)
    
    # Color code the cells based on the 'names' column
    colors <- colorRampPalette(c("blue", "red", "green", "yellow"))(length(unique(df$name)))
    color_mapping <- setNames(colors, unique(df$name))
    cell_style <- createStyle(fgFill = color_mapping[df_group$name])
    addStyle(wb, paste("Group", group), style = cell_style, rows = 2:(nrow(correlation_matrix)+1), cols = 2:(ncol(correlation_matrix)+1), gridExpand = TRUE)
  }
  
  # Save the workbook to an Excel file
  saveWorkbook(wb, paste0(category, "_correlations.xlsx"), overwrite = TRUE)
}
  
# Call the function for each category
create_excel(df, "pc_clusters")
create_excel(df, "labels")
create_excel(df, "clustering")


library(corrplot)
library(RColorBrewer)

# Function to create correlation plots for each group within a category
# Adjusted create_corrplot function
create_corrplot <- function(df, category) {
  # Get unique groups within the category
  groups <- unique(df[[category]])
  
  # For each group, compute the correlation matrix and create a correlation plot
  for (group in groups) {
    # Subset the data for this group
    df_group <- df %>% filter((!!sym(category)) == group)
    
    # Remove columns with zero standard deviation
    df_group <- df_group %>% 
      select(Tcell:DC, signature_1_dormant, residuals_dormant_ucellscore, all_of(jose_genes)) %>% 
      select_if(function(x) sd(x, na.rm=TRUE) != 0)
    
    # Check if the data frame is not empty
    if (nrow(df_group) > 0 && ncol(df_group) > 0) {
      # Compute the correlation matrix
      correlation_matrix <- df_group %>% 
        cor(use = "pairwise.complete.obs")
      
      # Create a correlation plot
      pdf(file = paste0(category, "_Group_", group, "_correlation.pdf"), width = 10, height = 10)
      corrplot(correlation_matrix, method = "color", col = colorRampPalette(c("blue", "white", "red"))(200), 
               title = paste0(category, "_Group_", group, "_correlation"), mar = c(0,0,1,0))
      dev.off()
    } else {
      message(paste("Skipping", category, "group", group, "- data frame is empty after removing columns with zero standard deviation."))
    }
  }
}


# Call the function for each category
create_corrplot(df, "pc_clusters")
create_corrplot(df, "labels")
create_corrplot(df, "clustering")



###########CORRELATION PLOTS AND EXCELL######### FIN


###########DAVID values AND EXCELL#########
# Group the data by 'labels', 'name', 'clustering', and 'pc_clusters', then calculate median values
df_grouped <- df %>%
  select(name, labels, clustering, pc_clusters, Tcell:DC, signature_1_dormant, all_of(jose_genes)) %>%
  group_by(name, labels, clustering, pc_clusters) %>%
  summarise(across(Tcell:DC, median, na.rm = TRUE),
            across(signature_1_dormant, median, na.rm = TRUE),
            across(all_of(jose_genes), median, na.rm = TRUE),
            .groups = "keep")

# Save 'pc_clusters' as last column
df_grouped <- df_grouped %>% relocate(pc_clusters, .after = last_col())

# Create two dataframes based on 'clustering' values
df_group1 <- df_grouped %>% filter(clustering %in% c(1,2,3,4))
df_group2 <- df_grouped %>% filter(clustering %in% c(5,6,7))

# Save dataframes to csv files
write.csv(df_group1, "./results/correlation/pc_groups/group1_dataframe.csv", row.names = FALSE)
write.csv(df_group2, "./results/correlation/pc_groups/group2_dataframe.csv", row.names = FALSE)

#####group them
# Create a new column 'sample_group' combining 'name' and 'labels'
df_group1$sample_group <- paste(df_group1$name, df_group1$labels, sep = "_")
df_group2$sample_group <- paste(df_group2$name, df_group2$labels, sep = "_")

# Group the data by 'sample_group' and calculate median values
df_group1_grouped <- df_group1 %>%
  group_by(sample_group) %>%
  summarise(across(Tcell:DC, median, na.rm = TRUE),
            across(signature_1_dormant, median, na.rm = TRUE),
            across(all_of(jose_genes), median, na.rm = TRUE),
            .groups = "keep")

df_group2_grouped <- df_group2 %>%
  group_by(sample_group) %>%
  summarise(across(Tcell:DC, median, na.rm = TRUE),
            across(signature_1_dormant, median, na.rm = TRUE),
            across(all_of(jose_genes), median, na.rm = TRUE),
            .groups = "keep")

# Save dataframes to csv files
write.csv(df_group1_grouped, "./results/correlation/pc_groups/group1_dataframe_grouped.csv", row.names = FALSE)
write.csv(df_group2_grouped, "./results/correlation/pc_groups/group2_dataframe_grouped.csv", row.names = FALSE)

##cor
# Install required packages
if (!require(corrplot)) install.packages("corrplot")
if (!require(gplots)) install.packages("gplots")

# Load required packages
library(corrplot)
library(gplots)

# Determine the columns with non-zero standard deviation
cols_to_keep1 <- apply(df_group1_grouped[,3:ncol(df_group1_grouped)], 2, sd) != 0
cols_to_keep2 <- apply(df_group2_grouped[,3:ncol(df_group2_grouped)], 2, sd) != 0

# Subset the data frames
df_grouped1_grouped <- df_group1_grouped[, c(TRUE, TRUE, cols_to_keep1)]
df_grouped2_grouped <- df_group2_grouped[, c(TRUE, TRUE, cols_to_keep2)]


# Compute correlations for df_grouped1_grouped and df_grouped2_grouped
corr1 <- cor(df_grouped1_grouped[,3:ncol(df_grouped1_grouped)], use = "pairwise.complete.obs")
corr2 <- cor(df_grouped2_grouped[,3:ncol(df_grouped2_grouped)], use = "pairwise.complete.obs")

hc_order <- hclust(dist(corr1))$order

# Set graphical parameters to make the labels smaller
par(cex=0.6)

# Generate correlation plot for df_grouped1_grouped
pdf("./results/correlation/pc_groups/Correlation_df_group1_grouped.pdf")
corrplot(corr1, method="color", col= colorRampPalette(c("blue", "white", "red"))(200),
         type="upper", addCoef.col = "black", tl.col="black", tl.srt=45, number.cex=0.70)
dev.off()

# Generate correlation plot for df_grouped2_grouped
pdf("./results/correlation/pc_groups/Correlation_df_group2_grouped.pdf")
corrplot(corr2, method="color", col= colorRampPalette(c("blue", "white", "red"))(200),
         type="upper", addCoef.col = "black", tl.col="black", tl.srt=45, number.cex=0.70)
dev.off()
##cor fin

##hierachical
# Perform hierarchical clustering on df_group1_grouped
dist_matrix <- dist(df_group1_grouped[,3:ncol(df_group1_grouped)])
hc <- hclust(dist_matrix)

# Plot the dendrogram
pdf("./results/correlation/pc_groups/Dendrogram_group1_grouped.pdf")
plot(hc, labels = df_group1_grouped$sample_group, main = "Hierarchical Clustering Dendrogram - Group1")
dev.off()

# Repeat for df_group2_grouped
dist_matrix <- dist(df_group2_grouped[,3:ncol(df_group2_grouped)])
hc <- hclust(dist_matrix)

# Plot the dendrogram
pdf("./results/correlation/pc_groups/Dendrogram_group2_grouped.pdf")
plot(hc, labels = df_group2_grouped$sample_group, main = "Hierarchical Clustering Dendrogram - Group2")
dev.off()
##hierachical fin

##cor rows
# Hierarchical clustering and Heatmap for df_group1_grouped
heatmap_data <- df_group1_grouped[,3:ncol(df_group1_grouped)]  # select only numeric columns
row_names <- df_group1_grouped$sample_group  # save row names
row.names(heatmap_data) <- row_names  # assign row names to heatmap data

# Create a heatmap
pdf("./results/correlation/pc_groups/Heatmap_group1_grouped.pdf")
heatmap.2(as.matrix(heatmap_data), 
          trace="none", 
          main="Heatmap - Group1",
          labRow=row_names,
          hclustfun = function(x) hclust(x, method="ward.D2"),
          distfun = function(x) dist(x, method="euclidean"))
dev.off()

# Repeat for df_group2_grouped

heatmap_data <- df_group2_grouped[,3:ncol(df_group2_grouped)]
row_names <- df_group2_grouped$sample_group
row.names(heatmap_data) <- row_names

pdf("./results/correlation/pc_groups/Heatmap_group2_grouped.pdf")
heatmap.2(as.matrix(heatmap_data), 
          trace="none", 
          main="Heatmap - Group2",
          labRow=row_names,
          hclustfun = function(x) hclust(x, method="ward.D2"),
          distfun = function(x) dist(x, method="euclidean"))
dev.off()

##cor_rows fin

##dendogram with pc_clusters
# Select only numeric columns and set row names as pc_clusters
heatmap_data <- df_group1[,3:ncol(df_group1)]
row.names(heatmap_data) <- df_group1$pc_clusters

# Compute a hierarchical clustering
hc <- hclust(dist(heatmap_data), method="ward.D2")

# Create a pdf for the dendrogram
pdf("./results/correlation/pc_groups/Dendrogram_group1.pdf")
plot(hc, main = "Hierarchical Clustering Dendrogram - df_grouped",
     xlab = "pc_clusters",
     sub = "",
     cex = 0.6)
dev.off()

heatmap_data <- df_group2[,3:ncol(df_group2)]
row.names(heatmap_data) <- df_group2$pc_clusters

# Compute a hierarchical clustering
hc <- hclust(dist(heatmap_data), method="ward.D2")

# Create a pdf for the dendrogram
pdf("./results/correlation/pc_groups/Dendrogram_group2.pdf")
plot(hc, main = "Hierarchical Clustering Dendrogram - df_grouped",
     xlab = "pc_clusters",
     sub = "",
     cex = 0.6)
dev.off()

##dendogram with pc_clusters fin

###########DAVID values AND EXCELL######### FIN







