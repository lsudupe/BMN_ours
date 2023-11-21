## SCRIPT: Human sample MM cleanning and QC

## 21.11.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(RColorBrewer)
color <- rev(brewer.pal(11,"Spectral"))

## Data
mm <- readRDS("./objects/sp/human/human_combined.rds")


## Spatial plots
for (i in 1:length(mm)){
  a <- mm[[i]]
  b <- names(mm[i])
  
  # count plot
  count <- "nCount_RNA"
  c <- c(min(a@meta.data[count]), max(a@meta.data[count]))
  p <- FeatureOverlay(a, features = count, ncols = 1, pt.size = 1.7, slot = "count",
                 value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c,
                         labels = c,
                         limits = c)
  # Save the plot to a PDF file
  pdf(paste0("./results/human/QC/mm/", b, "_ncount.pdf", sep = ""))
  print(p)
  dev.off()
  
  # feature plot
  feature <- "nFeature_RNA"
  c <- c(min(a@meta.data[feature]), max(a@meta.data[feature]))
  p <- FeatureOverlay(a, features = feature, ncols = 1, pt.size = 1.7, slot = "count",
                 value.scale = "all" ,cols = color) +
    scale_fill_gradientn(colours = color,
                         breaks = c,
                         labels = c,
                         limits = c)
  # Save the plot to a PDF file
  pdf(paste0("./results/human/QC/mm/", b, "_nfeature.pdf", sep = ""))
  print(p)
  dev.off()
  
}




