## SCRIPT: Clustering deconvolution results BM project

## 11.12.22 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(corrplot)
library(fpc)
library(dplyr)
library(clustertend)
library(factoextra)
library(gridExtra)
library(dendextend)

#Data---------------------------------
femur_M1 <- readRDS("./objects/card/heterogeneity/M1_fem_1C_subgroup_deco.rds")
femur_M3 <- readRDS("./objects/card/heterogeneity/M3_fem_1C_subgroup_deco.rds")

femur <- merge(femur_M1, y = c(femur_M3),  project = "BM")
femur

#Analysis--------------------------------
set.seed(20000)
meta <- femur@meta.data
types <- meta[,11:22]

matrix <- data.matrix(types, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/femur",filename = "both_cor.pdf"))
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

###hierarchycal plot
pdf(file.path("./results/endogram/femur",filename = "both_cor.pdf"))
get_clust_tendency(types, 2, graph=TRUE, gradient=list(low="red", mid="white", high="blue"))
dev.off()

###Optimal clustering

a <- fviz_nbclust(types, FUNcluster = kmeans, method = "silhouette") + theme_classic() 
b <- fviz_nbclust(types, FUNcluster = cluster::pam, method = "silhouette") + theme_classic() 
c <- fviz_nbclust(types, FUNcluster = cluster::clara, method = "silhouette") + theme_classic() 
d <- fviz_nbclust(types, FUNcluster = hcut, method = "silhouette") + theme_classic() 
e <- fviz_nbclust(types, FUNcluster = cluster::fanny, method = "silhouette") + theme_classic() 

pdf(file.path("./results/endogram/femur",filename = "both_minckuster.pdf"))
print(grid.arrange(a, b, c, d, e, ncol=2))
dev.off()

###hierarchical clustering
hc3 <- eclust(types, k=4, FUNcluster="hclust", hc_metric="euclidean", hc_method = "ward.D2")
hc3 %>% 
as.dendrogram()  -> dend


###ggplot automatic color scale
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 8
cols = gg_color_hue(n)
###

colors <- c("#F8766D", "#C49A00" ,"#53B400" ,"#00C094" ,"#00B6EB" ,"#A58AFF" ,"#FB61D7")
names(colors) = c("0", "1", "2", "3", "4","5", "6")
colors

col = ifelse(x$a == "Africa", "yellow",
             ifelse(x$a %in% c("North America", "South America"), "blue",
                    ifelse(x$a == "Asia", "green",
                           ifelse(x$a == "Europe", "lightblue",
                                  ifelse(x$a == "Oceania", "purple", "black")))))

col_1 = ifelse(meta$seurat_clusters == "0", "#F8766D",
             ifelse(meta$seurat_clusters == "1", "#C49A00",
                    ifelse(meta$seurat_clusters == "2", "#53B400", 
                           ifelse(meta$seurat_clusters == "3", "#00C094",
                                  ifelse(meta$seurat_clusters == "4", "#00B6EB",
                                         ifelse(meta$seurat_clusters == "5", "#A58AFF", 
                                                ifelse(meta$seurat_clusters == "6", "#FB61D7",
                                                       ifelse(meta$seurat_clusters == "7", "#FF61CC", 'white'))))))))
col = ifelse(meta$orig.ident, "grey", "gold")

col <- cbind(col_1, col_2)

pdf(file.path("./results/endogram/femur",filename = "both_dend.pdf"))
# Make the dendrogram
par(mar = c(10,2,1,1))
dend %>%
  set("labels_col", value = c("skyblue", "orange", "grey", "green"), k=4) %>%set("labels_cex", 0.1) %>%
  set("branches_k_color", value = c("skyblue", "orange", "grey", "green"), k = 4) %>%
  set("leaves_pch", 0.1)  %>% 
  plot()
colored_bars(colors = col, dend = dend, rowLabels = "seurat",sort_by_labels_order = FALSE)
dev.off()

# number of observations per cluster
hc_stats3 <- cluster.stats(types, hc3$cluster)
hc_stats3$cluster.size 

###ad clusters to spatial data
a <- as.factor(hc3[["cluster"]])
femur@meta.data[["clustering"]] <- a

###plot
b <- SetIdent(femur, value = femur@meta.data[["clustering"]])
pdf(file.path("./results/endogram/femur",filename = "both_spatial_hierarchical.pdf"))
print(SpatialDimPlot(b, combine = TRUE,label.size = 1.5, label = T, crop = TRUE, pt.size.factor = 10))
dev.off()

b <- SetIdent(femur, value = femur@meta.data[["seurat_clusters"]])
pdf(file.path("./results/endogram/femur",filename = "both_spatial_seurat.pdf"))
print(SpatialDimPlot(b, combine = TRUE,label.size = 1.5, label = T, crop = TRUE, pt.size.factor = 10))
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


##cor
matrix <- data.matrix(clust_3, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/femur",filename = "both_cluster3_cor.pdf"))
print(corrplot(M,
      type = 'upper',
               tl.col = "black",number.cex = 0.75))
dev.off()


####Extract hierarchycal clustering porcentages

#######proportions loop
value <- as.vector(unique(types_b$clustering))
lista <- list()

for (i in value){
#select cluster of interest rows
value_1 <- types_b[grepl(i, types_b$clustering),]
value_1$clustering <- NULL

###create list to add content
proportions <- c()
  for (o in colnames(value_1)){ 
    proportions <- c(proportions, (sum(value_1[[o]])*100/nrow(value_1)))
  }
name <- paste('cluster:',i,sep='')
lista[[name]] <- proportions

}



columns = colnames(types_b)
df = data.frame(matrix(nrow = length(value), ncol = length(columns))) 
colnames(df) = columns
rownames(df) = value
df

for (i in colnames(df)){ 
  df
  proportions <- c(proportions, (sum(value_1[[o]])*100/nrow(value_1)))
}
name <- paste('cluster:',i,sep='')
lista[[name]] <- proportions

}

## dataframe
pro_df_ <- data.frame(BM_proportions)
rownames(pro_df_) <- z[1:12]
write.csv(pro_df_, file = paste0("./results/CARD/heterogeneity/",v,"pro.csv"))



