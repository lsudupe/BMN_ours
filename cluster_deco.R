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

#Data--------------------------------
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
hc3 <- eclust(types, k=6, FUNcluster="hclust", hc_metric="euclidean", hc_method = "ward.D2")
hc3 %>% 
as.dendrogram()  -> dend

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
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

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

#######CORRELATION 
types_b

##cor cluster6
clust_6 <- types_b[grepl("6", types_b$clustering),]
clust_6$clustering <- NULL
matrix <- data.matrix(clust_6, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/st/",filename = "cor_cluster6_square.pdf"))
print(corrplot(M, method = 'square', title="Cluster 6 cell type correlation",
               col=colorRampPalette(c("blue","white","red"))(100),
               order = 'original', type = 'lower', diag = FALSE,
               mar=c(0,0,1,0)))
dev.off()


##cor cluster5
clust_5 <- types_b[grepl("5", types_b$clustering),]
clust_5$clustering <- NULL
matrix <- data.matrix(clust_5, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/st/",filename = "cor_cluster5_square.pdf"))
print(corrplot(M, method = 'square', title="Cluster 5 cell type correlation",
               col=colorRampPalette(c("blue","white","red"))(100),
               order = 'original', type = 'lower',diag = FALSE,
               mar=c(0,0,1,0)))
dev.off()

##cor cluster4
clust_4 <- types_b[grepl("4", types_b$clustering),]
clust_4$clustering <- NULL
matrix <- data.matrix(clust_4, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/st/",filename = "cor_cluster4_square.pdf"))
print(corrplot(M, method = 'square', title="Cluster 4 cell type correlation",
               col=colorRampPalette(c("blue","white","red"))(100),
               order = 'original', type = 'lower', diag = FALSE,
               mar=c(0,0,1,0)))
dev.off()

types_b$clustering <- NULL

matrix <- data.matrix(types_b, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/st/",filename = "cor_al_square.pdf"))
print(corrplot(M, method = 'square', title="MM samples cell type correlation",
               col=colorRampPalette(c("blue","white","red"))(100),
               order = 'original', type = 'lower', diag = FALSE,
               mar=c(0,0,1,0)))
dev.off()

#######CORRELATION fin

####Extract hierarchycal clustering porcentages
#######proportions loop
value <- as.vector(unique(types_b$clustering))
#value <- as.vector(unique(types_b$sample))

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

write.csv(df, "./celltypes.csv", row.names=TRUE)

# Basic piechart
df["celltype"] <- as.vector(rownames(df))
pdf(file.path("./results/endogram/st",filename = "clustering_percentages_piechart.pdf"))
ggplot(df, aes(x="", y=df$`sample:M1_fem_1C`, fill=df$celltype)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = df$`sample:M1_fem_1C`),
            position = position_stack(vjust = 0.5)) +
  coord_polar("y", start=0)
dev.off()


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
                 cluster6 = c(df$`cluster:6`)
                 )
DFtall <- DF %>% gather(key = Cluster, value = Value, cluster1:cluster6)
DFtall

pdf(file.path("./results/endogram/st",filename = "clustering_percentages.pdf"))
ggplot(DFtall, aes(Cluster, Value, fill = group)) + geom_col(position = "dodge")
dev.off()

# Stacked + percent
#http://rstudio-pubs-static.s3.amazonaws.com/5312_98fc1aba2d5740dd849a5ab797cc2c8d.html
pdf(file.path("./results/endogram/st",filename = "clustering_percentages_barplot.pdf"))
ggplot(DFtall, aes(fill=group, y=Value, x=Cluster)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

ggplot(DFtall, aes(fill=group, y=Value, x=Cluster)) + 
  geom_bar(colour = "black", position = "stack") + scale_fill_brewer(palette = "RdBu") + 
  labs(title = "BrBG")

pdf(file.path("./results/endogram/st",filename = "clustering_percentages_barplot_prueba_2.pdf"))
ggplot(DFtall, aes(fill=group, y=Value, x=Cluster)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_brewer(palette = "PuOr")
dev.off()


##borrar no healthy
new <- DFtall[!grepl("cluster5", DFtall$Cluster),]
new <- new[!grepl("cluster6", new$Cluster),]

# Stacked + percent
pdf(file.path("./results/endogram/st",filename = "clustering_percentages_barplot_onlyhealthy.pdf"))
ggplot(new, aes(fill=group, y=Value, x=Cluster)) + 
  geom_bar(position="fill", stat="identity")
dev.off()


