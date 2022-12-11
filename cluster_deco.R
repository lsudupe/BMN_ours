## SCRIPT: Clustering deconvolution results BM project

## 11.12.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(corrplot)
library(fpc)

#Data---------------------------------
femur_M1 <- readRDS("./objects/card/heterogeneity/M1_fem_1C_subgroup_deco.rds")
femur_M3 <- readRDS("./objects/card/heterogeneity/M3_fem_1C_subgroup_deco.rds")


#Analysis--------------------------------
set.seed(20000)
meta <- femur@meta.data
types <- meta[,11:22]

matrix <- data.matrix(types, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/femur",filename = "M1_cor.pdf"))
print(corrplot(M, method = "number", number.cex = 0.75, order="hclust"))
dev.off()

vars <- c(colnames(types))
Outliers <- c()
for(i in vars){
  max <- quantile(types[,i], 0.75) + (IQR(types[,i]) * 1.5 )
  min <- quantile(types[,i], 0.25) - (IQR(types[,i]) * 1.5 )
  idx <- which(types[,i] < min | types[,i] > max)
  print(paste(i, length(idx), sep=' ')) # printing variable and number of potential outliers 
  Outliers <- c(Outliers, idx) 
}

Outliers
##plot te outliers
par(mfrow=c(2,2))
colnames <- colnames(types[,c(2:5,9:11)])
for (i in colnames) {
  plot(types[,i], main = paste("Plot of ", i), ylab = i)
}

###Prediagnostic
library(clustertend)
library(factoextra)
get_clust_tendency(types, 2, graph=TRUE, gradient=list(low="red", mid="white", high="blue"))

###Optimal clustering
library(gridExtra)
a <- fviz_nbclust(types, FUNcluster = kmeans, method = "silhouette") + theme_classic() 
b <- fviz_nbclust(types, FUNcluster = cluster::pam, method = "silhouette") + theme_classic() 
c <- fviz_nbclust(types, FUNcluster = cluster::clara, method = "silhouette") + theme_classic() 
d <- fviz_nbclust(types, FUNcluster = hcut, method = "silhouette") + theme_classic() 
e <- fviz_nbclust(types, FUNcluster = cluster::fanny, method = "silhouette") + theme_classic() 
grid.arrange(a, b, c, d, e, ncol=2)


###hierarchical clustering
hc1 <- eclust(types, k=3, FUNcluster="hclust", hc_metric="euclidean", hc_method = "complete")
plot(hc1, cex=0.6, hang=-1, main = "Dendrogram of HAC")
rect.hclust(hc1, k=3, border='red')

hc3 <- eclust(types, k=5, FUNcluster="hclust", hc_metric="euclidean", hc_method = "ward.D2")
plot(hc3, cex=0.6, hang=-1, main = "Dendrogram of HAC")
rect.hclust(hc3, k=3, border='red')
# number of observations per cluster
hc_stats3 <- cluster.stats(types, hc3$cluster)
hc_stats3$cluster.size 

###ad clusters to spatial data
a <- as.factor(hc3[["cluster"]])
femur_M1@meta.data[["clustering"]] <- a

###plot
b <- SetIdent(femur_M1, value = femur_M1@meta.data[["clustering"]])
pdf(file.path("./results/endogram/femur",filename = "M1_spatial_hierarchical.pdf"))
print(SpatialDimPlot(b, combine = FALSE,label.size = 1.5, label = T, crop = TRUE, pt.size.factor = 10))
dev.off()

b <- SetIdent(femur_M1, value = femur_M1@meta.data[["seurat_clusters"]])
pdf(file.path("./results/endogram/femur",filename = "M1_spatial_seurat.pdf"))
print(SpatialDimPlot(b, combine = FALSE,label.size = 1.5, label = T, crop = TRUE, pt.size.factor = 10))
dev.off()

###subset the data in hierarchical clustering
meta_b <- b@meta.data
types_b <- meta_b[,11:23]

clust_1 <- types_b[grepl("1", types_b[,13]),]
clust_1 <- clust_1[,1:12]
clust_2 <- types_b[grepl("2", types_b[,13]),]
clust_2 <- clust_2[,1:12]
clust_3 <- types_b[grepl("3", types_b[,13]),]
clust_3 <- clust_3[,1:12]
clust_4 <- types_b[grepl("4", types_b[,13]),]
clust_4 <- clust_4[,1:12]
clust_5 <- types_b[grepl("5", types_b[,13]),]
clust_5 <- clust_5[,1:12]

##cor
matrix <- data.matrix(clust_5, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/femur",filename = "cluster5_cor.pdf"))
print(corrplot(M, method = "number", number.cex = 0.75, order="hclust"))
dev.off()
