## SCRIPT: Tcells genes individual plot  BM project

## 11.07.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(dplyr)
library(UCell)
library(tidyr)
library(RColorBrewer)
library(tidyverse)

#Data---------------------------------
all <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

# Extract the gene list from the Seurat object
seurat_genes <- rownames(all@assays$RNA@counts)

#pc_markers
tcells16 <- c("Cd8a", "Cd4", "Foxp3", "Il2ra", "Pdcd1", "Lag3", "Tigit", "Ifng", 
                   "Gzma", "Gzmb", "Gzmk", "Ctla4", "Cxcr3", "Tnfrsf9", "Icos", "Tnfrsf4")

tcells8 <- c("Cd4", "Foxp3", "Tigit", "Ctla4", "Cxcr3", "Tnfrsf9", "Icos", "Tnfrsf4")

#CD45,CD56,CD117
# Check which genes_to_check are in seurat_genes
genes_present <- tcells8 %in% seurat_genes


###score
## Add UCellScore
#16genes
vector<- ScoreSignatures_UCell(all@assays[["RNA"]]@counts, features = list(tcells16))
all@meta.data[["signature_1_tcells"]] <- as.vector(vector)

color <- brewer.pal(11,"Spectral")
color <- rev(color)

p <- FeatureOverlay(all, features = c("signature_1_tcells"),sampleids = 1:6, ncols = 2, pt.size = 0.7, 
                    value.scale = "all" ,cols = color)

pdf(paste("./results/tcells/mouse/genes16_all.pdf",sep=""))
print(p)
dev.off()

#8genes
vector<- ScoreSignatures_UCell(all@assays[["RNA"]]@counts, features = list(tcells8))
all@meta.data[["signature_1_tcells8"]] <- as.vector(vector)

color <- brewer.pal(11,"Spectral")
color <- rev(color)

p <- FeatureOverlay(all, features = c("signature_1_tcells8"),sampleids = 1:6, ncols = 2, pt.size = 0.7, 
                    value.scale = "all" ,cols = color)

pdf(paste("./results/tcells/mouse/genes8_all.pdf",sep=""))
print(p)
dev.off()

##individual genes
for (i in tcells8){
  a <- i
  p <- FeatureOverlay(all, features = i,sampleids = 1:6, ncols = 2, pt.size = 0.7, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/tcells/mouse/", i,"_tcell.pdf",sep=""))
  print(p)
  dev.off()
  
  plot <- VlnPlot(object = all, features = i, split.by = "clustering")
  
  pdf(paste("./results/tcells/mouse/", i,"_tcell_violin.pdf",sep=""), width = 10, height = 7)
  print(plot)
  dev.off()
  
  Idents(object = all) <- "clustering"
  plot <- VlnPlot(object = all, features = i)
  
  pdf(paste("./results/tcells/mouse/", i,"_tcell_violin_ALL.pdf",sep=""), width = 10, height = 7)
  print(plot)
  dev.off()
  
  #boxplot
  # Fetch data for the current feature
  data <- FetchData(all, vars = i)
  
  # Add cluster identities to the data
  data$cluster <- all@meta.data[["clustering"]]
  
  # Convert to long format data frame for ggplot
  data_long <- data %>% pivot_longer(cols = starts_with(i),
                                     names_to = "Feature",
                                     values_to = "Expression") %>% 
    filter(Expression != 0)
  
  # Plot boxplot using ggplot2
  p <- ggplot(data_long, aes(x = cluster, y = Expression, fill = cluster)) + 
    geom_boxplot() +
    labs(title = i) +
    theme_minimal() +
    xlab("Cluster") +
    ylab("Expression Level")

  pdf(paste("./results/tcells/mouse/", i,"_tcell_boxplot.pdf",sep=""), width = 10, height = 7)
  print(p)
  dev.off()
  
}

Idents(object = all) <- "clustering"
VlnPlot(object = all, features = c("Cd4"))

