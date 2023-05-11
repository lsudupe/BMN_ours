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

meta <- se@meta.data
types <- meta[,12:18]

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
types_b <- meta_b[,12:18]
types_b["clustering"] <- as.vector(b@meta.data[["clustering"]])

###violin plot
ggplot(types_b, aes(x = clustering, y = MM_MIC)) +
  geom_violin(scale = "width", trim = FALSE, show.legend = FALSE) +
  stat_summary(fun.data = mean_sdl, color = "red", geom = "pointrange", position = position_dodge(0.9)) +
  stat_summary(fun = median, color = "blue", geom = "point", position = position_dodge(0.9)) +
  geom_boxplot(width = 0.1, outlier.color = "black", outlier.shape = 16, outlier.size = 1) +
  theme_minimal() +
  labs(x = "Clustering", y = "Cell Type", title = "Violin Plot of Plasma by Clustering")

##distribution plot
pdf(file.path("./results/endogram/st",filename = "clustering_distribution.pdf"))
types_b %>%
  pivot_longer(cols = -clustering) %>%
  mutate(Cluster = factor(clustering)) %>%
  ggplot(aes(x = name, y = value, fill = Cluster)) + geom_violin()
dev.off()

#######CORRELATION 
types_b

##cor cluster7
clust_7 <- types_b[grepl("7", types_b$clustering),]
clust_7$clustering <- NULL
matrix <- data.matrix(clust_7, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/st/",filename = "cor_cluster7_square.pdf"))
print(corrplot(M, method = 'square', title="Cluster 7 cell type correlation",
               col=colorRampPalette(c("blue","white","red"))(100),
               order = 'original', type = 'lower', diag = FALSE,
               mar=c(0,0,1,0)))
dev.off()

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

##cor cluster3
clust_3 <- types_b[grepl("3", types_b$clustering),]
clust_3$clustering <- NULL
matrix <- data.matrix(clust_3, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/endogram/st/",filename = "cor_cluster3_square.pdf"))
print(corrplot(M, method = 'square', title="Cluster 3 cell type correlation",
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
types_b["clustering"] <- as.vector(b@meta.data[["clustering"]])
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

# Stacked + percent
#http://rstudio-pubs-static.s3.amazonaws.com/5312_98fc1aba2d5740dd849a5ab797cc2c8d.html
pdf(file.path("./results/endogram/st",filename = "clustering_percentages_barplot.pdf"))
ggplot(DFtall, aes(fill=group, y=Value, x=Cluster)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

library(RColorBrewer)
cell_type_colors <- brewer.pal(length(rows), "RdBu")
cells_order <- c("Bcell", "DC", "Erythroblasts","MM_MIC", "Monocytes","Neutrophils","Tcell")
cell_type_color_map <- setNames(cell_type_colors, cells_order)

pdf(file.path("./results/endogram/st",filename = "clustering_percentages_barplot_colors.pdf"))
ggplot(DFtall, aes(fill=group, y=Value, x=Cluster)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = cell_type_color_map)
dev.off()


##borrar no healthy
new <- DFtall[!grepl("cluster5", DFtall$Cluster),]
new <- new[!grepl("cluster6", new$Cluster),]
new <- new[!grepl("cluster7", new$Cluster),]
new <- new[!grepl("cluster3", new$Cluster),]


# Stacked + percent
pdf(file.path("./results/endogram/st",filename = "clustering_percentages_barplot_onlyhealthy.pdf"))
ggplot(new, aes(fill=group, y=Value, x=Cluster)) + 
  geom_bar(position="fill", stat="identity")
dev.off()

#########Each cluster plots
p1 <- ST.FeaturePlot(se, features = "clustering", indices = 1, 
                     split.labels = T, pt.size = 0.3) & theme(plot.title = element_blank(), strip.text = element_blank())
p2 <- ST.FeaturePlot(se, features = "clustering", indices = 2, 
                     split.labels = T, pt.size = 0.3) & theme(plot.title = element_blank(), strip.text = element_blank())

p3 <- ST.FeaturePlot(se, features = "clustering", indices = 5, 
                     split.labels = T, pt.size = 0.3) & theme(plot.title = element_blank(), strip.text = element_blank())

p4 <- ST.FeaturePlot(se, features = "clustering", indices = 6, 
                     split.labels = T, pt.size = 0.3) & theme(plot.title = element_blank(), strip.text = element_blank())

p5 <- ST.FeaturePlot(se, features = "clustering", indices = 3, 
                     split.labels = T, pt.size = 0.3) & theme(plot.title = element_blank(), strip.text = element_blank())

pdf(file.path("./results/endogram/st",filename = "each_cluster.pdf"))
cowplot::plot_grid(p1, p2, p3, p4,p5, ncol = 5)
dev.off()

#########Each cluster plots FIN

#########DE analysis with pseudobulk approach
##DA.DB fb
#lm <- lm(meta$ratio_stand ~ meta$geneSetFB, data =meta)
#residuals <- lm$residuals
#a@meta.data[["residualsDADB_after_FB"]] <- residuals
library(edgeR)
library(sva)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(11,"Spectral")
color <- rev(color)

###1.subset data, only cluster 4 and 6
de <- SetIdent(se, value = se@meta.data[["clustering"]])
de.subset <- SubsetSTData(de, idents = c("4", "5","6"))
de.subset45 <- SubsetSTData(de, idents = c("4","5"))

FeatureOverlay(de, features = "clustering", sampleids = 5:6, pt.size = 1.3)
FeatureOverlay(de.subset, features = "clustering", sampleids = 3:4, pt.size = 1.3)
FeatureOverlay(de.subset45, features = "clustering", sampleids = 5:6, pt.size = 1.3)


###2.create pseudobulk
cts <- AggregateExpression(de.subset, 
                           group.by = c("clustering", "name"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$RNA
df <- as.data.frame(cts)
df$`4_M3_F_1C` <- NULL
df$`4_M3_fem_1C` <- NULL
df$`4_M1_fem_1C`<- NULL
df$`4_M9_F2_1C` <- NULL
colnames(df)<- c("cluster4_M2_F_2B","cluster4_M8_F2_1C","cluster6_M2_F_2B","cluster6_M8_F2_1C")
cts <- as.matrix(df)

###3. pseudobulk analysis

# Named vector of sample names
sample <- c("M2_F_2B","M8_F2_1C")
sample = as.factor(c(rep(sample,2)))
clusters <- as.factor(c("cluster4","cluster4","cluster6","cluster6"))
mycols = data.frame(row.names = colnames(df),sample, clusters)
mycols
# Our area of interest
mycols$clusters <- relevel(mycols$clusters, "cluster4")
# Create dsd object
dsd = DESeqDataSetFromMatrix(countData = cts, colData = mycols, design = ~ clusters) 

# Run DESeq2 differential expression analysis
dds <- DESeq(dsd)
dds$clusters <- relevel(dds$clusters, "cluster4")
dataset <- counts(dds,  normalized = TRUE)

## batch correction
adjusted_counts <- ComBat_seq(dataset, batch=mycols$sample,
                              group=mycols$clusters)
# Normalize the counts
normalize <- NormalizeData(adjusted_counts)

## EDGER
y <- DGEList(counts=adjusted_counts,group=mycols$clusters)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
designGO <- model.matrix(~ 0+ mycols$clusters)

# contrast of interest
colnames(designGO)<-c("cluster4","cluster6")
y <- estimateDisp(y,designGO)
con <- makeContrasts(cluster4 - cluster6, levels=designGO)
fit <- glmQLFit(y,designGO)
qlf <- glmQLFTest(fit,contrast=con)

# cluster4 - cluster6
qlf.4vs6 <- glmQLFTest(fit,contrast=con[,"cluster4 - cluster6"])
qlf.4vs6_matrix <- as.data.frame(topTags(qlf.4vs6,n=10700,sort.by = "PValue"))
genes_of_interest_qlf.4vs6 <- row.names(qlf.4vs6_matrix[qlf.4vs6_matrix[,"FDR"]<0.006 & qlf.4vs6_matrix$logFC>1,])
top.genes  <- unique(genes_of_interest_qlf.4vs6)

## prepare data
normalize_matrix <- as.matrix(normalize)
matrix <- normalize_matrix[top.genes, ]

matrix_genes <- data.frame(normalize_matrix) %>%
  rownames_to_column(var = "genes") %>%
  dplyr::filter(genes %in% top.genes)


matrix_genes_ <- data.frame(matrix_genes[,-1], row.names = matrix_genes[,1])

# Run pheatmap using the metadata data frame for the annotation
pdf("./results/DE/st/cluster4vscluster6.pdf")
print(pheatmap(matrix_genes_, 
               #color = heat_colors, 
               cluster_rows = T, 
               cluster_cols = F,
               show_rownames = T,
               cutree_cols = 4,
               annotation = mycols[, c("clusters", "sample")],
               scale = "row", 
               fontsize_row = 5, 
               height = 20)) 
dev.off()

#########DE analysis with pseudobulk approach FIN


#########clusterprofile START
#ClusterProfiler
library("clusterProfiler")
library("org.Mm.eg.db")
library("AnnotationHub")


markers <- FindAllMarkers(de.subset, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

##################################################################
#Subsetting top 100 markers with adjusted p values lower than .05#
##################################################################
top100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)

df <- top100pval[,7:6]
dfsample <- split(df$gene,df$cluster)
length(dfsample)

pdf("./results/DE/st/cluster456_de.pdf")
DoHeatmap(de.subset, features = top100pval$gene,disp.min = -2, disp.max = 2)
dev.off()

#The output of length(dfsample) returns how many clusters you have
#Here there at 9 clusters (0, 1, 2, 3, 4, 5, 6, 7 and 8)
#I'm sure there's a better way but you have to make a line like below for each cluster

dfsample$`4` = bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`5` = bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`6` = bitr(dfsample$`6`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

#do the same here, a line like below for each cluster
genelist <- list("4" = dfsample$`4`$ENTREZID,
                 "5" = dfsample$`5`$ENTREZID,
                 "6" = dfsample$`6`$ENTREZID)

GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
ck <- setReadable(GOclusterplot, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

head(ck) 
pdf("./results/DE/st/go.456.pdf")
print(dotplot(GOclusterplot))
dev.off()

pdf("./results/DE/st/cnetplot.456.pdf")
cnetplot(ck)
dev.off()

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
dotplot(KEGGclusterplot)

#########clusterprofile FIN
