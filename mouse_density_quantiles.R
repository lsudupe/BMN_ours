## SCRIPT: Mouse quantile

## 21.11.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(UCell)
library(STutility)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(openxlsx)
color <- rev(brewer.pal(11,"Spectral"))

## Functions
source("./regress_out_function.R")

## Data
mouse <- readRDS("./objects/sp/se_hierarchical_signatures.rds")
Idents(object = mouse) <- "name"
mouse <- SubsetSTData(mouse, idents = c("M1_fem_1C", "M2_F_2B", "M8_F2_1C", "M9_F2_1C"))

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
#mouse <- regress_out(mouse, "MM_MIC")
mouse <- SCTransform(mouse, vars.to.regress = "MM_MIC", new.assay.name ="SCT_PC")

## DE

Idents(mouse) <- "Quantile"
# Differential expression for Q3 vs Q4
de_Q3_Q4 <- FindMarkers(mouse, ident.1 = "Q3", ident.2 = "Q4", min.pct = 0.25)
# Differential expression for Q2 vs Q3
de_Q2_Q3 <- FindMarkers(mouse, ident.1 = "Q2", ident.2 = "Q3", min.pct = 0.25)

## Comparation plots
## Volcan
prepare_for_volcano <- function(de_results, top_n = 50) {
  de_results$logP <- -log10(de_results$p_val_adj)
  de_results$logFC <- de_results$avg_log2FC
  de_results$significant <- de_results$p_val_adj < 0.05 & abs(de_results$logFC) > log2(1.5)
  
  top_genes_to_label <- de_results[de_results$significant, ]
  top_genes_to_label <- top_genes_to_label[order(top_genes_to_label$logP, decreasing = TRUE), ][1:top_n, ]
  
  list(data = de_results, labels = top_genes_to_label)
}

data_Q3_Q4 <- prepare_for_volcano(de_Q3_Q4)
data_Q2_Q3 <- prepare_for_volcano(de_Q2_Q3)


create_volcano_plot <- function(data, title) {
  ggplot(data$data, aes(x=logFC, y=logP)) +
    geom_point(aes(color=significant), alpha=0.5) +
    geom_text_repel(data = data$labels, aes(label=rownames(data$labels)), size = 3, box.padding = 0.35, point.padding = 0.5) +
    theme_minimal() +
    scale_color_manual(values=c("grey", "red")) +
    labs(title=title, x="Log2 Fold Change", y="-Log10 Adjusted p-value") +
    geom_vline(xintercept=0, linetype="dashed", color = "black", alpha=0.7) +
    theme(legend.position="none")
}

# Create and plot for "Q3 vs Q4"
volcano_Q3_Q4 <- create_volcano_plot(data_Q3_Q4, "Volcano Plot: Q3 vs Q4")
volcano_Q3_Q4
# Create and plot for "Q2 vs Q3"
volcano_Q2_Q3 <- create_volcano_plot(data_Q2_Q3, "Volcano Plot: Q2 vs Q3")
volcano_Q2_Q3

## Violin
VlnPlot(mouse, features = c("Cldn7", "Slamf9"))


## Save excell
# Create a new workbook
wb <- createWorkbook()

# Add sheets with the DE results
de_Q3_Q4$Gene <- rownames(de_Q3_Q4)
de_Q2_Q3$Gene <- rownames(de_Q2_Q3)

addWorksheet(wb, "DE_Q3_Q4")
writeData(wb, sheet = "DE_Q3_Q4", de_Q3_Q4)

addWorksheet(wb, "DE_Q2_Q3")
writeData(wb, sheet = "DE_Q2_Q3", de_Q2_Q3)

# Save the workbook to a file
saveWorkbook(wb, "./DE_quantiles.xlsx", overwrite = TRUE)

