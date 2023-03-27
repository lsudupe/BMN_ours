## SCRIPT: QC of the ST spatial objects BMN project

## 27.03.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

## Read data
se <- readRDS("./objects/sc/integrated/se_deco.rds")



## Plots

library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(11,"Spectral")
color <- rev(color)
  
  
p <- ST.FeaturePlot(se, features = "nFeature_RNA", ncol = 2 ,pt.size = 0.65, show.sb = FALSE, palette = "Spectral")
p <-FeatureOverlay(se, features = c("nFeature_RNA"), sampleids = 1:6,pt.size = 0.40,ncol = 2 , 
                   value.scale = "all" ,cols = color)



pdf(file.path("./results/ST/QC/",filename = "nfeature.pdf"))
print(p & theme_bw())
dev.off()
