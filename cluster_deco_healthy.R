## SCRIPT: Clustering deconvolution results BM project healthy samples

## 28.03.23 Laura Sudupe , git @lsudupe

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
se <- readRDS("./objects/sc/integrated/se_deco_healthy.rds")


#Analysis--------------------------------
set.seed(20000)

meta <- se@meta.data
types <- meta[,12:20]

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
pdf(file.path("./results/endogram/st/healthy/",filename = "both_cor_.pdf"))
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
hc3 <- eclust(types, k=4, FUNcluster="hclust", hc_metric="euclidean", hc_method = "ward.D2")
hc3 %>% 
as.dendrogram()  -> dend

# number of observations per cluster
hc_stats3 <- cluster.stats(types, hc3$cluster)
hc_stats3$cluster.size 

###ad clusters to spatial data
a <- as.factor(hc3[["cluster"]])
#femur@meta.data[["clustering"]] <- a
se@meta.data[["clustering"]] <- a


##save object
saveRDS(se, "./objects/heterogeneity/healthy/se_hierarchical.rds")
se <- readRDS("./objects/heterogeneity/healthy/se_hierarchical.rds")

FeatureOverlay(se, features = c("Ighm"), ncol = 2,pt.size = 0.7)
###plot
b <- SetIdent(se, value = se@meta.data[["clustering"]])
pdf(file.path("./results/endogram/st/healthy/",filename = "both_spatial_hierarchical.pdf"))
print(FeatureOverlay(b, features = "clustering", sampleids = 1:2, pt.size = 1.3))
dev.off()

b <- SetIdent(se, value = se@meta.data[["seurat_clusters"]])
pdf.options(encoding = )
pdf(file.path("./results/endogram/st",filename = "both_spatial_seurat.pdf"))
print(FeatureOverlay(b, features = "seurat_clusters", sampleids = 1:6, ncols = 2,pt.size = 0.7))
dev.off()

###subset the data in hierarchical clustering
meta_b <- b@meta.data
types_b <- meta_b[,12:20]
types_b["clustering"] <- as.vector(b@meta.data[["clustering"]])
types_b["sample"] <- as.vector(b@meta.data[["name"]])


##distribution plot
pdf(file.path("./results/endogram/st/healthy/",filename = "clustering_distribution.pdf"))
types_b %>%
  pivot_longer(cols = -clustering) %>%
  mutate(Cluster = factor(clustering)) %>%
  ggplot(aes(x = name, y = value, fill = Cluster)) + geom_violin()
dev.off()

types_only <- types_b
types_only$clustering <- NULL

matrix <- data.matrix(types_only, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/st/healthy/",filename = "healthy_cor.pdf"))
print(corrplot(M, method = 'square', 
               col=colorRampPalette(c("lightblue","white","red"))(100),
               order = 'FPC', type = 'lower', diag = FALSE))
dev.off()


##cor BY CLUSTER
matrix <- data.matrix(clust_6, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/st/",filename = "cluster6_cor.pdf"))
print(corrplot(M, method = 'square', 
               col=colorRampPalette(c("lightblue","white","red"))(100),
               order = 'FPC', type = 'lower', diag = FALSE))
dev.off()

##cor BY CLUSTER FIN

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
df = data.frame(matrix(nrow = length(colnames(value_1)), ncol = 0)) 
rownames(df) = colnames(value_1)

##add list values to df
for (i in 1:length(lista)){
  a <- lista[[i]]
  df[, ncol(df) + 1] <- a
  names(df)[ncol(df)] <- names(lista[i])
}

write.csv(df, "./celltypes.csv", row.names=TRUE)

#df["celltype"] <- as.vector(rownames(df))

# Basic piechart
pdf(file.path("./results/endogram/st/healthy/",filename = "clustering_percentages_piechart_cluster1.pdf"))
ggplot(df, aes(x="", y=df$`cluster:1`, fill=df$celltype)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label = df$`cluster:1`),
            position = position_stack(vjust = 0.5)) +
  coord_polar("y", start=0)
dev.off()


# Stacked + percent
pdf(file.path("./results/endogram/st/healthy/",filename = "clustering_percentages_barplot_onlyhealthy.pdf"))
ggplot(df, aes(fill=group, y=Value, x=Cluster)) + 
  geom_bar(position="fill", stat="identity")
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
                 cluster4 = c(df$`cluster:4`)
                 )
DFtall <- DF %>% gather(key = Cluster, value = Value, cluster1:cluster4)
DFtall

pdf(file.path("./results/endogram/st/healthy/",filename = "clustering_percentages.pdf"))
ggplot(DFtall, aes(Cluster, Value, fill = group)) + geom_col(position = "dodge")
dev.off()

# Stacked + percent
pdf(file.path("./results/endogram/st/healthy/",filename = "clustering_percentages_barplot_onlyhealthy.pdf"))
ggplot(DFtall, aes(fill=group, y=Value, x=Cluster)) + 
  geom_bar(position="fill", stat="identity")
dev.off()


######VIOLIN PLOT 

se <- readRDS("./objects/heterogeneity/healthy/se_hierarchical.rds")
b <- SetIdent(se, value = se@meta.data[["clustering"]])
###subset the data in hierarchical clustering
meta_b <- b@meta.data
types_b <- meta_b[,12:20]
types_b["sample"] <- as.vector(b@meta.data[["name"]])

##distribution plot
pdf(file.path("./results/endogram/st/healthy/",filename = "sample_distribution.pdf"))
types_b %>%
  pivot_longer(cols = -sample) %>%
  mutate(Sample = factor(sample)) %>%
  ggplot(aes(x = name, y = value, fill = Sample)) + geom_violin()
dev.off()

value <- as.vector(unique(types_b$sample))

lista <- list()
for (i in value){
  #select cluster of interest rows
  value_1 <- types_b[grepl(i, types_b$sample),]
  value_1$sample <- NULL
  
  ###create list to add content
  proportions <- c()
  for (o in colnames(value_1)){ 
    proportions <- c(proportions, (sum(value_1[[o]])*100/nrow(value_1)))
  }
  name <- paste('sample:',i,sep='')
  lista[[name]] <- proportions
}


######create a df 
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
                 M3_F_1C = c(df$`sample:M3_F_1C`),
                 M3_fem_1C = c(df$`sample:M3_fem_1C`)
)
DFtall <- DF %>% gather(key = sample, value = Value, M3_F_1C:M3_fem_1C)
DFtall


# Stacked + percent
pdf(file.path("./results/endogram/st/healthy/",filename = "sample_percentages_barplot.pdf"))
ggplot(DFtall, aes(fill=group, y=Value, x=sample)) + 
  geom_bar(position="fill", stat="identity")
dev.off()



######VIOLIN PLOT FIN


