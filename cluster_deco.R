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
library(tidyr)
library(STutility)

#Data---------------------------------
femur_M1 <- readRDS("./objects/card/heterogeneity/M1_fem_1C_subgroup_deco.rds")
femur_M3 <- readRDS("./objects/card/heterogeneity/M3_fem_1C_subgroup_deco.rds")

femur <- merge(femur_M1, y = c(femur_M3),  project = "BM")
femur

se <- readRDS("./objects/sc/integrated/se_deco.rds")


#Analysis--------------------------------
set.seed(20000)
meta <- femur@meta.data
types <- meta[,10:19]

meta <- se@meta.data
types <- meta[,12:21]

matrix <- data.matrix(types, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/st",filename = "both_cor.pdf"))
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
colnames <- colnames(types[,c(1:10)])
for (i in colnames) {
  plot(types[,i], main = paste("Plot of ", i), ylab = i)
}

###Prediagnostic

###hierarchycal plot
pdf(file.path("./results/endogram/st",filename = "both_cor_.pdf"))
get_clust_tendency(types, 2, graph=TRUE, gradient=list(low="red", mid="white", high="blue"))
dev.off()

###Optimal clustering

a <- fviz_nbclust(types, FUNcluster = kmeans, method = "silhouette") + theme_classic() 
b <- fviz_nbclust(types, FUNcluster = cluster::pam, method = "silhouette") + theme_classic() 
c <- fviz_nbclust(types, FUNcluster = cluster::clara, method = "silhouette") + theme_classic() 
d <- fviz_nbclust(types, FUNcluster = hcut, method = "silhouette") + theme_classic() 
e <- fviz_nbclust(types, FUNcluster = cluster::fanny, method = "silhouette") + theme_classic() 

pdf(file.path("./results/endogram/st",filename = "both_minckuster.pdf"))
print(grid.arrange(a, b, c, d, ncol=2))
dev.off()

###hierarchical clustering
hc3 <- eclust(types, k=7, FUNcluster="hclust", hc_metric="euclidean", hc_method = "ward.D2")
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

#col <- cbind(col_1, col_2)

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
#femur@meta.data[["clustering"]] <- a
se@meta.data[["clustering"]] <- a


##save object
saveRDS(se, "./objects/heterogeneity/se_hierarchical.rds")

###plot
b <- SetIdent(se, value = se@meta.data[["clustering"]])
pdf(file.path("./results/endogram/st",filename = "both_spatial_hierarchical.pdf"))
print(FeatureOverlay(b, features = "clustering", sampleids = 1:6, ncols = 2,pt.size = 0.7))
dev.off()

b <- SetIdent(se, value = se@meta.data[["seurat_clusters"]])
pdf(file.path("./results/endogram/st",filename = "both_spatial_seurat.pdf"))
print(FeatureOverlay(b, features = "seurat_clusters", sampleids = 1:6, ncols = 2,pt.size = 0.7))
dev.off()

###subset the data in hierarchical clustering
meta_b <- b@meta.data
types_b <- meta_b[,12:21]
types_b["clustering"] <- as.vector(b@meta.data[["clustering"]])

##distribution plot
pdf(file.path("./results/endogram/st",filename = "clustering_distribution.pdf"))
types_b %>%
  pivot_longer(cols = -clustering) %>%
  mutate(Cluster = factor(clustering)) %>%
  ggplot(aes(x = name, y = value, fill = Cluster)) + geom_violin()
dev.off()

clust_1 <- types_b[grepl("1", types_b$clustering),]
clust_1$clustering <- NULL
clust_2 <- types_b[grepl("2", types_b$clustering),]
clust_2$clustering <- NULL
clust_3 <- types_b[grepl("3", types_b$clustering),]
clust_3$clustering <- NULL
clust_4 <- types_b[grepl("4", types_b$clustering),]
clust_4$clustering <- NULL


##cor
matrix <- data.matrix(clust_4, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/femur",filename = "both_cluster4_cor.pdf"))
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

######create a df 
rows = colnames(value_1)
df = data.frame(matrix(nrow = length(colnames(value_1)), ncol = 0)) 
rownames(df) = colnames(value_1)

##add list values to df
for (i in 1:length(lista)){
  a <- lista[[i]]
  df[, ncol(df) + 1] <- a
  names(df)[ncol(df)] <- names(lista[i])
}

#df <- t(df)
df <- cbind(celltype = rownames(df), df)
rownames(df) <- 1:nrow(df)

library(tidyr)
library(ggplot2)
DF <- data.frame(group = c(df$celltype),
                 cluster1 = c(df$`cluster:1`),
                 cluster2 = c(df$`cluster:2`),
                 cluster3 = c(df$`cluster:3`),
                 cluster4 = c(df$`cluster:4`),
                 cluster5 = c(df$`cluster:5`),
                 cluster6 = c(df$`cluster:6`),
                 cluster7 = c(df$`cluster:7`)
                 )
DFtall <- DF %>% gather(key = Cluster, value = Value, cluster1:cluster7)
DFtall

pdf(file.path("./results/endogram/st",filename = "clustering_percentages.pdf"))
ggplot(DFtall, aes(Cluster, Value, fill = group)) + geom_col(position = "dodge")
dev.off()



