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
library(hrbrthemes)

#Data---------------------------------
all <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

# Extract the gene list from the Seurat object
seurat_genes <- rownames(all@assays$RNA@counts)

#pc_markers
tcells16 <- c("Cd8a", "Cd4", "Foxp3", "Il2ra", "Pdcd1", "Lag3", "Tigit", "Ifng", 
                   "Gzma", "Gzmb", "Gzmk", "Ctla4", "Cxcr3", "Tnfrsf9", "Icos", "Tnfrsf4")

tcells8 <- c("Cd4", "Foxp3", "Tigit", "Ctla4", "Cxcr3", "Tnfrsf9", "Icos", "Tnfrsf4")

candidate <- c("Pdcd1", "Ctla4", "Tnfrsf9", "Havcr2", "Tox", "Tigit", "Wars", "Rsad2",
               "Mcm7", "Mx1", "Ndfip2", "Enosf1", "Ccdc141", "Stmn1", "Ttn", "Faslg",
               "Mcm5", "Nab1", "Phlda1", "Mcm3", "Pcna", "Gapdh", "Oasl", "Ifi44l",
               "Tbc1d4", "Slc43a3", "Pam", "Ccl3", "Acp5", "Oas3", "Cd38", "Tnfs10",
               "Gbp2", "Kif20b", "Ctsb")

#CD45,CD56,CD117
# Check which genes_to_check are in seurat_genes
genes_present <- tcells8 %in% seurat_genes


###score
## Add UCellScore
#candidate
vector<- ScoreSignatures_UCell(all@assays[["RNA"]]@counts, features = list(candidate))
all@meta.data[["signature_1_candidate"]] <- as.vector(vector)

p <- FeatureOverlay(all, features = c("signature_1_candidate"),sampleids = 1:6, ncols = 2, pt.size = 0.7, 
                    value.scale = "all" ,cols = color)

pdf(paste("./results/tcells/mouse/candidate_all.pdf",sep=""))
print(p)
dev.off()

Idents(object = all) <- "clustering"
plot <- VlnPlot(object = all, features = "signature_1_candidate")

pdf(paste("./results/tcells/mouse/candidate_violin.pdf",sep=""), width = 10, height = 7)
print(plot)
dev.off()

##Regressout Tcell value
## Regress out FB values
meta <- all@meta.data
lm <- lm(meta$signature_1_candidate ~ meta$Tcell, data =meta)
residuals <- lm$residuals
all@meta.data[["residuals_exhausted"]] <- residuals

p <- FeatureOverlay(all, features = c("residuals_exhausted"),sampleids = 1:6, ncols = 2, pt.size = 0.7, 
                    value.scale = "all" ,cols = color)

pdf(paste("./results/tcells/mouse/residuals_exhausted_all.pdf",sep=""))
print(p)
dev.off()

Idents(object = all) <- "clustering"
plot <- VlnPlot(object = all, features = "residuals_exhausted")

pdf(paste("./results/tcells/mouse/residuals_exhausted_violin.pdf",sep=""), width = 10, height = 7)
print(plot)
dev.off()

plot <- VlnPlot(all, "residuals_exhausted", group.by = "name", split.by = "clustering") +
  geom_boxplot()

pdf(paste("./results/tcells/mouse/residuals_exhausted_violin_bysample.pdf",sep=""), width = 10, height = 7)
print(plot)
dev.off()

plot <- VlnPlot(all, "residuals_exhausted", group.by = "clustering", split.by = "clustering") +
  geom_boxplot()

pdf(paste("./results/tcells/mouse/residuals_exhausted_violin_allbox.pdf",sep=""), width = 10, height = 7)
print(plot)
dev.off()


###Tcells
plot <- VlnPlot(all, "Tcell", group.by = "name", split.by = "clustering") +
  geom_boxplot()

pdf(paste("./results/tcells/mouse/Tcell_violin_bysample.pdf",sep=""), width = 10, height = 7)
print(plot)
dev.off()

###scater plot
df <- all@meta.data
ggplot(df, aes(x=signature_1_candidate, y=MM_MIC, color=clustering)) + 
  geom_point(size=0.5) +
  geom_smooth(method=lm , color="red") +
  theme_ipsum()

plot <- VlnPlot(all, "signature_1_candidate", group.by = "clustering", split.by = "clustering") +
  geom_boxplot()

pdf(paste("./results/tcells/mouse/candidate_violin.pdf",sep=""), width = 10, height = 7)
print(plot)
dev.off()

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

specific_genes <- c("Lag3", "Tigit", "Cd69", "Cd83", "Btla", "Slamf6", 
                    "Ctla4")

#Cd43, Pd1, Tim3 no

for (i in specific_genes){
  a <- i
  p <- FeatureOverlay(all, features = i,sampleids = 1:6, ncols = 2, pt.size = 0.7, 
                      value.scale = "all" ,cols = color)
  
  pdf(paste("./results/tcells/mouse/specific_genes/", i,"_tcell.pdf",sep=""))
  print(p)
  dev.off()
  
  plot <- VlnPlot(object = all, features = i, split.by = "clustering")
  
  pdf(paste("./results/tcells/mouse//specific_genes/", i,"_tcell_violin.pdf",sep=""), width = 10, height = 7)
  print(plot)
  dev.off()
  
  Idents(object = all) <- "clustering"
  plot <- VlnPlot(object = all, features = i)
  
  pdf(paste("./results/tcells/mouse/specific_genes/", i,"_tcell_violin_ALL.pdf",sep=""), width = 10, height = 7)
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
  
  pdf(paste("./results/tcells/mouse/specific_genes/", i,"_tcell_boxplot.pdf",sep=""), width = 10, height = 7)
  print(p)
  dev.off()
  
}

