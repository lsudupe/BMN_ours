## SCRIPT: Plot gradients results BM project healthy samples

## 02.04.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(dplyr)
library(hrbrthemes)
library(extrafont)
library(ggthemes)
library(viridis)
library("scales")

#Data---------------------------------
se <- readRDS("./objects/heterogeneity/healthy/se_hierarchical.rds")
x <- se

##replace all values below 10% with 0
meta <- x@meta.data
a <- meta[,12:20]
a[a < 0.1] <- 0
b <- apply(a, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

meta[,12:20] <- b
x@meta.data <- meta

#x.subset <- SubsetSTData(x, expression = Tcell >= 0.1)
#se.subset <- SubsetSTData(se, expression = nFeature_RNA >= 2000)

####cor
matrix <- data.matrix(b, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/ST/gradient/healthy/prueba_cor.pdf"))
print(corrplot(M, method = 'square', title="Cell type correlation",
               col=colorRampPalette(c("blue","white","red"))(100),
               order = 'original', type = 'lower', diag = FALSE,
               mar=c(0,0,1,0)))
dev.off()

####cor fin

####spatial plots


pdf(file.path("./results/ST/gradient/healthy/prueba.pdf"))
print(ST.FeaturePlot(x, features = c("Erythroblasts", "Neutrophils"),sampleids = 1, pt.size = 0.7,ncol = 1 , 
                     value.scale = "all" , palette = "Spectral"))
dev.off()

p1 <- ST.FeaturePlot(x, features = c("Erythroblasts", "Neutrophils"),sampleids = 1, pt.size = 0.5,ncol = 1 , 
                     value.scale = "all" , palette = "Spectral")
p2 <- VlnPlot(x, features = c("Erythroblasts", "Neutrophils"), ncol = 2, group.by = "clustering")


pdf(file.path("./results/ST/gradient/healthy/prueba_2.pdf"))
print(p1 + p2)
dev.off()
