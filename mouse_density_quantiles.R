## SCRIPT: Mouse quantile

## 21.11.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(UCell)
library(STutility)
library(dplyr)
library(RColorBrewer)
library("viridis")
color <- rev(brewer.pal(11,"Spectral"))

## Functions
source("./regress_out_function.R")

## Data
mouse <- readRDS("./objects/sp/se_hierarchical_signatures.rds")

## Quantiles
metadata <- mouse@meta.data

# Calculate quartiles
Q1 <- quantile(metadata$MM_MIC, 0.25)
Q2 <- quantile(metadata$MM_MIC, 0.5)
Q3 <- quantile(metadata$MM_MIC, 0.75)
Q4 <- quantile(metadata$MM_MIC, 1)

# Function to assign quartile
assignQuartile <- function(value) {
  if (value <= Q1) {
    return('Q1')
  } else if (value <= Q2) {
    return('Q2')
  } else if (value <= Q3) {
    return('Q3')
  } else {
    return('Q4')
  }
}

# Apply the function to create a new column 'Quantile'
metadata$Quantile <- sapply(metadata$MM_MIC, assignQuartile)
# Update the Seurat object
mouse@meta.data <- metadata

## Quantile spatial plot
FeatureOverlay(mouse, features = "Quantile", ncols = 2, pt.size = 1,
               value.scale = "all" ,cols = color, sampleids = 1:6) 


## Regress out pc
mouse <- regress_out(mouse, "MM_MIC")

## DE

Idents(mouse) <- "Quantile"
all_markers <- FindAllMarkers(mouse)

top_markers_per_cluster <- all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.005,avg_log2FC > 2) %>%
  slice_head(n = 5) 

## Comparation plots
# Create a heatmap for the top genes
DoHeatmap(mouse, features = top_markers_per_cluster$gene, group.by = c("Quantile", "clustering")) 


VlnPlot(mouse, features = c("Sdc1", "Sel1l", "Ly6d", "Pon3", "Cp"))

VlnPlot(mouse, features = c("Igha", "Eef1a1", "B2m", "Sec11c", "Hist1h1d"))

## Volcano
# Select results for a specific quantile, e.g., Q1
quantile_results <- all_markers[all_markers$cluster == "Q4", ]

# Prepare the data for plotting
quantile_results$logP <- -log10(quantile_results$p_val_adj) # Using adjusted p-values
quantile_results$logFC <- quantile_results$avg_log2FC

# Adding a column to label genes with significant changes
threshold_pval <- 0.05
threshold_logFC <- log2(1.5)

quantile_results$significant <- quantile_results$p_val_adj < threshold_pval & 
  abs(quantile_results$logFC) > threshold_logFC

top_genes_to_label <- quantile_results[quantile_results$significant, ]
top_genes_to_label <- top_genes_to_label[order(top_genes_to_label$logP, decreasing = TRUE), ][1:10, ] # Top 10 as an example

# Creating the volcano plot with labels
library(ggrepel)

ggplot(quantile_results, aes(x=logFC, y=logP)) +
  geom_point(aes(color=significant), alpha=0.5) +
  geom_text_repel(data = top_genes_to_label, aes(label=gene), size = 3, box.padding = 0.35, point.padding = 0.5) +
  theme_minimal() +
  scale_color_manual(values=c("grey", "red")) +
  labs(title="Volcano Plot: Q4",
       x="Log2 Fold Change",
       y="-Log10 Adjusted p-value") +
  geom_vline(xintercept=0, linetype="dashed", color = "black", alpha=0.7) +
  theme(legend.position="none")
