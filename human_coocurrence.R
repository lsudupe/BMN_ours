## SCRIPT: Human sample pc, tex, reads co-ocurrence

## 11.10.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(Signac)
library(RColorBrewer)
library(STutility)
library(UCell)
library(readr)
library(dplyr)
library(gridExtra)
library(tibble)

save_dir <- "/Users/medinils/Desktop/BMN_ours/results/human/coocurrence/"

#pc_markers
pc_mm_markers <- c("SHH", "DHH", "IHH", "PTCH1", "PTCH2", "SMO", "SUFU", "GLI1", 
                   "GLI2", "GLI3", "CD19", "CD44", "CXCR4", "KLF4", "CD28", "CD33", "CD27")

human_tcell_exhausted <- c("PDCD1", "CTLA4", "TNFRSF9", "HAVCR2", "TOX", "TIGIT", "WARS", "RSAD2",
                           "MCM7", "MX1", "NDFIP2", "ENOSF1", "CCDC141", "STMN1", "TTN", "FASLG",
                           "MCM5", "NAB1", "PHLDA1", "MCM3", "PCNA", "GAPDH", "OASL", "IFI44L",
                           "TBC1D4", "SLC43A3", "PAM", "CCL3", "ACP5", "OAS3", "CD38", "TNFSF10",
                           "GBP2", "KIF20B", "CTSB")

human <- readRDS("./objects/sp/human/human_combined.rds")


# split human
B08041 <- human[["BM_human_AP-B08041_"]]
B08805 <- human[["BM_human_AP-B08805"]]
B10395 <- human[["BM_B10395"]]

lista <- c(B08041, B08805, B10395)
names(lista) <- c("B08041","B08805", "B10395")

###Enrichment score
color <- brewer.pal(11,"Spectral")
color <- rev(color)

post <- c()
for (i in 1:length(lista)){
  a <- lista[[i]]
  b <- names(lista[i])
  
  # Add UCellScore
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(pc_mm_markers))
  a@meta.data[["signature_1_pc"]] <- as.vector(vector)
    
  vector<- ScoreSignatures_UCell(a@assays[["RNA"]]@counts, features = list(human_tcell_exhausted))
  a@meta.data[["signature_1_tex"]] <- as.vector(vector)
  
  # Create the plot
  p <- FeatureOverlay(a, features = c("nCount_RNA"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste(save_dir, b,"_nCount_RNA.pdf",sep=""))
  print(p)
  dev.off()
  
  p <- FeatureOverlay(a, features = c("nFeature_RNA"), ncols = 1, pt.size = 1.1, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste(save_dir, b,"_nFeature_RNA.pdf",sep=""))
  print(p)
  dev.off()
  
  
  
  post[[length(post) + 1]] <- a
}

names(post) <- c("B08041","B08805", "B10395")


#######

# Extract metadata from each Seurat object
meta_B08041 <- post$B08041@meta.data
meta_B08805 <- post$B08805@meta.data
meta_B10395 <- post$B10395@meta.data

# Calculate and print correlations
correlations <- list()

correlations$B08041 <- list(
  signature_1_pc_nCount_RNA = cor(meta_B08041$signature_1_pc, meta_B08041$nCount_RNA),
  signature_1_tex_nCount_RNA = cor(meta_B08041$signature_1_tex, meta_B08041$nCount_RNA),
  signature_1_pc_nFeature_RNA = cor(meta_B08041$signature_1_pc, meta_B08041$nFeature_RNA),
  signature_1_tex_nFeature_RNA = cor(meta_B08041$signature_1_tex, meta_B08041$nFeature_RNA)
)

correlations$B08805 <- list(
  signature_1_pc_nCount_RNA = cor(meta_B08805$signature_1_pc, meta_B08805$nCount_RNA),
  signature_1_tex_nCount_RNA = cor(meta_B08805$signature_1_tex, meta_B08805$nCount_RNA),
  signature_1_pc_nFeature_RNA = cor(meta_B08805$signature_1_pc, meta_B08805$nFeature_RNA),
  signature_1_tex_nFeature_RNA = cor(meta_B08805$signature_1_tex, meta_B08805$nFeature_RNA)
)

correlations$B10395 <- list(
  signature_1_pc_nCount_RNA = cor(meta_B10395$signature_1_pc, meta_B10395$nCount_RNA),
  signature_1_tex_nCount_RNA = cor(meta_B10395$signature_1_tex, meta_B10395$nCount_RNA),
  signature_1_pc_nFeature_RNA = cor(meta_B10395$signature_1_pc, meta_B10395$nFeature_RNA),
  signature_1_tex_nFeature_RNA = cor(meta_B10395$signature_1_tex, meta_B10395$nFeature_RNA)
)

print(correlations)

# Scatter plots for visualization
plot_scatter <- function(data, x, y, title) {
  ggplot(data, aes_string(x = x, y = y)) + 
    geom_point(alpha = 0.5) +
    geom_smooth(method = 'lm', col = 'red') + # Linear regression line
    ggtitle(title) +
    theme_minimal()
}

# For B08041
plot1 <- plot_scatter(meta_B08041, "signature_1_pc", "nCount_RNA", "B08041: signature_1_pc vs nCount_RNA")
plot2 <- plot_scatter(meta_B08041, "signature_1_tex", "nCount_RNA", "B08041: signature_1_tex vs nCount_RNA")
plot3 <- plot_scatter(meta_B08041, "signature_1_pc", "nFeature_RNA", "B08041: signature_1_pc vs nFeature_RNA")
plot4 <- plot_scatter(meta_B08041, "signature_1_tex", "nFeature_RNA", "B08041: signature_1_tex vs nFeature_RNA")

# Function to save the scatter plots
save_plot <- function(plot, filename) {
  full_path <- paste0(save_dir, filename)
  ggsave(filename = full_path, plot = plot, width = 6, height = 4)
}

# Save the scatter plots for B08041
save_plot(plot1, "B08041_signature_1_pc_vs_nCount_RNA.pdf")
save_plot(plot2, "B08041_signature_1_tex_vs_nCount_RNA.pdf")
save_plot(plot3, "B08041_signature_1_pc_vs_nFeature_RNA.pdf")
save_plot(plot4, "B08041_signature_1_tex_vs_nFeature_RNA.pdf")


# For B08805
plot1 <- plot_scatter(meta_B08805, "signature_1_pc", "nCount_RNA", "B08041: signature_1_pc vs nCount_RNA")
plot2 <- plot_scatter(meta_B08805, "signature_1_tex", "nCount_RNA", "B08041: signature_1_tex vs nCount_RNA")
plot3 <- plot_scatter(meta_B08805, "signature_1_pc", "nFeature_RNA", "B08041: signature_1_pc vs nFeature_RNA")
plot4 <- plot_scatter(meta_B08805, "signature_1_tex", "nFeature_RNA", "B08041: signature_1_tex vs nFeature_RNA")


# Save the scatter plots for B08805
save_plot(plot1, "B08805_signature_1_pc_vs_nCount_RNA.pdf")
save_plot(plot2, "B08805_signature_1_tex_vs_nCount_RNA.pdf")
save_plot(plot3, "B08805_signature_1_pc_vs_nFeature_RNA.pdf")
save_plot(plot4, "B08805_signature_1_tex_vs_nFeature_RNA.pdf")

# For B10395
plot1 <- plot_scatter(meta_B10395, "signature_1_pc", "nCount_RNA", "B08041: signature_1_pc vs nCount_RNA")
plot2 <- plot_scatter(meta_B10395, "signature_1_tex", "nCount_RNA", "B08041: signature_1_tex vs nCount_RNA")
plot3 <- plot_scatter(meta_B10395, "signature_1_pc", "nFeature_RNA", "B08041: signature_1_pc vs nFeature_RNA")
plot4 <- plot_scatter(meta_B10395, "signature_1_tex", "nFeature_RNA", "B08041: signature_1_tex vs nFeature_RNA")


# Save the scatter plots for B08041
save_plot(plot1, "B10395_signature_1_pc_vs_nCount_RNA.pdf")
save_plot(plot2, "B10395_signature_1_tex_vs_nCount_RNA.pdf")
save_plot(plot3, "B10395_signature_1_pc_vs_nFeature_RNA.pdf")
save_plot(plot4, "B10395_signature_1_tex_vs_nFeature_RNA.pdf")

# Spearman's Rank Correlation function
spearman_correlation <- function(data, var1, var2) {
  cor.test(data[[var1]], data[[var2]], method = "spearman")
}

# Get Spearman's correlation for each pair in each sample
samples <- list(B08041 = meta_B08041, B08805 = meta_B08805, B10395 = meta_B10395)
cor_results <- list()

for (sample_name in names(samples)) {
  cor_results[[sample_name]] <- list(
    signature_1_pc_nCount_RNA = spearman_correlation(samples[[sample_name]], "signature_1_pc", "nCount_RNA"),
    signature_1_tex_nCount_RNA = spearman_correlation(samples[[sample_name]], "signature_1_tex", "nCount_RNA"),
    signature_1_pc_nFeature_RNA = spearman_correlation(samples[[sample_name]], "signature_1_pc", "nFeature_RNA"),
    signature_1_tex_nFeature_RNA = spearman_correlation(samples[[sample_name]], "signature_1_tex", "nFeature_RNA")
  )
}

# Print the Spearman's correlation results
print(cor_results)

# 2D density plot function
density_plot <- function(data, x, y, title) {
  ggplot(data, aes_string(x = x, y = y)) + 
    geom_density_2d_filled() +
    ggtitle(title) +
    theme_minimal()
}

# Save plots function
save_dir <- "/Users/medinils/Desktop/BMN_ours/results/human/coocurrence/"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

save_density_plot <- function(data, x, y, filename) {
  plot <- density_plot(data, x, y, filename)
  full_path <- paste0(save_dir, filename, ".pdf")
  ggsave(filename = full_path, plot = plot, width = 6, height = 4, device = "pdf")
}

# Save the density plots for all samples
for (sample_name in names(samples)) {
  save_density_plot(samples[[sample_name]], "signature_1_pc", "nCount_RNA", paste0(sample_name, "_signature_1_pc_vs_nCount_RNA_density"))
  save_density_plot(samples[[sample_name]], "signature_1_tex", "nCount_RNA", paste0(sample_name, "_signature_1_tex_vs_nCount_RNA_density"))
  save_density_plot(samples[[sample_name]], "signature_1_pc", "nFeature_RNA", paste0(sample_name, "_signature_1_pc_vs_nFeature_RNA_density"))
  save_density_plot(samples[[sample_name]], "signature_1_tex", "nFeature_RNA", paste0(sample_name, "_signature_1_tex_vs_nFeature_RNA_density"))
}

