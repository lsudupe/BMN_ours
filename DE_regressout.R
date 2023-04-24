## SCRIPT: DE results WITH REGRESS OUT BM project

## 11.04.23 Laura Sudupe , git @lsudupe

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
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
color <- brewer.pal(11,"Spectral")
color <- rev(color)

#Data--------------------------------
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")
b <- SetIdent(se, value = se@meta.data[["clustering"]])
DefaultAssay(b) <- "RNA"

####DE by samples and PC as a covariate
##divide by sample
Idents(object = b) <- "name"
name <- unique(b@meta.data[["name"]])
objects <- c()


for (i in name){
  a <- SubsetSTData(b, idents = i)
  ## add object to list
  objects[[length(objects) + 1]] <- a
  
}
names(objects) <- name
objects$se <- "se"
objects[["se"]] <- se

## separate list
#list2env(objects,envir=.GlobalEnv)

#set colors
n = 6
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(n)
cluster_color_map <- setNames(cols, unique(se@meta.data[["clustering"]]))


####Extract the matrix, regress out PC percentage, add the new assay to the objects
objects_new <- c()

for (i in 1:length(objects)){
  
  a <- objects[[i]]
  b <- names(objects[i])
  a <- ScaleData(a)
  ###subset the data 
  meta <- a@meta.data
  types <- meta[,12:21]
  types["clustering"] <- as.vector(a@meta.data[["clustering"]])
  
  ###extract counts
  matrix <- as.data.frame(a@assays[["RNA"]]@counts)
  
  matrix_t <- t(matrix)
  matrix_t <- as.data.frame(matrix_t)
  matrix_t["plasma_value"] <- types$MM_MIC
  matrix_final <- as.data.frame(t(matrix_t))
  
  ###REGRESS OUT
  data <- matrix_final
  # Get the number of genes and spots
  genes <- nrow(data) - 1 # Subtract 1 to exclude the last row with plasma cell percentage values
  spots <- ncol(data)
  # Convert the data frame to a matrix
  data_matrix <- as.matrix(data)
  # Initialize the residuals matrix with the same dimensions as the data matrix (without the last row)
  residuals_matrix <- matrix(0, nrow = genes, ncol = spots)
  # Loop through each gene
  for (i in 1:genes) {
    # Create a linear model to regress out the plasma cell percentage
    model <- lm(data_matrix[i, 1:spots] ~ data_matrix[genes + 1, 1:spots])
    # Save the residuals in the corresponding row of the residuals_matrix
    residuals_matrix[i, ] <- model$residuals
  }
  # Replace the original data_matrix with the residuals
  data_matrix[1:genes, 1:spots] <- residuals_matrix
  # Convert the matrix back to a data frame
  data <- as.data.frame(data_matrix)
  
  ###REGRESS OUT FIN
  
  ### add new infor to the object
  regress <- CreateAssayObject(data)
  a@assays[["regress"]] <- regress
  a@assays$regress@key <- "regress_"
  
  ### scale
  DefaultAssay(a) <- "regress"
  #a <- ScaleData(a)
  
  ## PLOTS
  ##Cd81
  DefaultAssay(a) <- "RNA"
  p <- VlnPlot(object = a, features = c("Cd81"), group.by = "clustering") + ggtitle("normal") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  DefaultAssay(a) <- "regress"
  pdf(paste("./results/ST/jose_angel/", b,"_Cd81_normal.pdf",sep=""))
  print(p)
  dev.off()
  p <-FeatureOverlay(se, features = c("Cd81"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                     value.scale = "all" ,cols = color)
  pdf(paste("./results/ST/jose_angel/", b,"_Cd81_spatial.pdf",sep=""))
  print(p)
  dev.off()
  
  DefaultAssay(a) <- "regress"
  p <- VlnPlot(object = a, features = c("Cd81"), group.by = "clustering") + ggtitle("regress out") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  pdf(paste("./results/ST/jose_angel/", b,"_Cd81_regress.pdf",sep=""))
  print(p)
  dev.off()

  ##Cd44
  DefaultAssay(a) <- "RNA"
  p <- VlnPlot(object = a, features = c("Cd44"), group.by = "clustering") + ggtitle("normal") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  DefaultAssay(a) <- "regress"
  pdf(paste("./results/ST/jose_angel/", b,"_Cd44_normal.pdf",sep=""))
  print(p)
  dev.off()
  p <-FeatureOverlay(se, features = c("Cd44"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                     value.scale = "all" ,cols = color)
  pdf(paste("./results/ST/jose_angel/", b,"_Cd44_spatial.pdf",sep=""))
  print(p)
  dev.off()
  
  DefaultAssay(a) <- "regress"
  p <- VlnPlot(object = a, features = c("Cd44"), group.by = "clustering") + ggtitle("regress out") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  pdf(paste("./results/ST/jose_angel/", b,"_Cd44_regress.pdf",sep=""))
  print(p)
  dev.off()
  
  ##Flna
  DefaultAssay(a) <- "RNA"
  p <- VlnPlot(object = a, features = c("Flna"), group.by = "clustering") + ggtitle("normal") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  DefaultAssay(a) <- "regress"
  pdf(paste("./results/ST/jose_angel/", b,"_Flna_normal.pdf",sep=""))
  print(p)
  dev.off()
  p <-FeatureOverlay(se, features = c("Flna"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                     value.scale = "all" ,cols = color)
  pdf(paste("./results/ST/jose_angel/", b,"_Flna_spatial.pdf",sep=""))
  print(p)
  dev.off()
  
  DefaultAssay(a) <- "regress"
  p <- VlnPlot(object = a, features = c("Flna"), group.by = "clustering") + ggtitle("regress out") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  pdf(paste("./results/ST/jose_angel/", b,"_Flna_regress.pdf",sep=""))
  print(p)
  dev.off()
  ##Mki67
  DefaultAssay(a) <- "RNA"
  p <- VlnPlot(object = a, features = c("Mki67"), group.by = "clustering") + ggtitle("normal") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  DefaultAssay(a) <- "regress"
  pdf(paste("./results/ST/jose_angel/", b,"_Mki67_normal.pdf",sep=""))
  print(p)
  dev.off()
  p <-FeatureOverlay(se, features = c("Mki67"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                     value.scale = "all" ,cols = color)
  pdf(paste("./results/ST/jose_angel/", b,"_Mki67_spatial.pdf",sep=""))
  print(p)
  dev.off()
  
  DefaultAssay(a) <- "regress"
  p <- VlnPlot(object = a, features = c("Mki67"), group.by = "clustering") + ggtitle("regress out") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  pdf(paste("./results/ST/jose_angel/", b,"_Mki67_regress.pdf",sep=""))
  print(p)
  dev.off()
  ##Pcna
  DefaultAssay(a) <- "RNA"
  p <- VlnPlot(object = a, features = c("Pcna"), group.by = "clustering") + ggtitle("normal") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  DefaultAssay(a) <- "regress"
  pdf(paste("./results/ST/jose_angel/", b,"_Pcna_normal.pdf",sep=""))
  print(p)
  dev.off()
  p <-FeatureOverlay(se, features = c("Pcna"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                     value.scale = "all" ,cols = color)
  pdf(paste("./results/ST/jose_angel/", b,"_Pcna_spatial.pdf",sep=""))
  print(p)
  dev.off()
  
  DefaultAssay(a) <- "regress"
  p <- VlnPlot(object = a, features = c("Pcna"), group.by = "clustering") + ggtitle("regress out") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  pdf(paste("./results/ST/jose_angel/", b,"_Pcna_regress.pdf",sep=""))
  print(p)
  dev.off()
  ##Xbp1
  DefaultAssay(a) <- "RNA"
  p <- VlnPlot(object = a, features = c("Xbp1"), group.by = "clustering") + ggtitle("normal") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  DefaultAssay(a) <- "regress"
  pdf(paste("./results/ST/jose_angel/", b,"_Xbp1_normal.pdf",sep=""))
  print(p)
  dev.off()
  p <-FeatureOverlay(se, features = c("Xbp1"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                     value.scale = "all" ,cols = color)
  pdf(paste("./results/ST/jose_angel/", b,"_Xbp1_spatial.pdf",sep=""))
  print(p)
  dev.off()
  
  DefaultAssay(a) <- "regress"
  p <- VlnPlot(object = a, features = c("Xbp1"), group.by = "clustering") + ggtitle("regress out") +geom_violin(trim=FALSE)+
    scale_fill_manual(values = cluster_color_map) +
    theme_minimal()
  pdf(paste("./results/ST/jose_angel/", b,"_Xbp1_regress.pdf",sep=""))
  print(p)
  dev.off()
  ## PLOTS FIN

  saveRDS(a,file = paste0("./objects/sp/regress_out/",b,"_regressout_ST.rds"))
  
  ## add object to list
  objects_new[[length(objects_new) + 1]] <- a
  
}
names(objects_new) <- name

saveRDS(objects_new,file = "./objects/sp/regress_out/list_regressout_ST.rds")

####integrated object DE
DefaultAssay(b) <- "SCT"
DefaultAssay(b) <- "RNA"
markers <- FindAllMarkers(b,min.pct = 0.1, logfc.threshold = 0.25)
cluster <- subset(markers, p_val_adj < 0.05 & 0.05 < avg_log2FC)
top100 <- cluster %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

pdf("./results/DE/st/se_all_clusters_doheatmap_rna.pdf")
DoHeatmap(b, assay = "RNA", features = top100$gene)#,disp.min = -2, disp.max = 2)
dev.off()

#cluster profiler
df <- top100[,7:6]
dfsample <- split(df$gene,df$cluster)
length(dfsample)

dfsample$`0` = bitr(dfsample$`0`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`1` = bitr(dfsample$`1`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`2` = bitr(dfsample$`2`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`3` = bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`4` = bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`5` = bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`6` = bitr(dfsample$`6`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

#do the same here, a line like below for each cluster
genelist <- list("0" = dfsample$`0`$ENTREZID,
                 "1" = dfsample$`1`$ENTREZID,
                 "2" = dfsample$`2`$ENTREZID,
                 "3" = dfsample$`3`$ENTREZID,
                 "4" = dfsample$`4`$ENTREZID,
                 "5" = dfsample$`5`$ENTREZID,
                 "6" = dfsample$`6`$ENTREZID)

GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
ck <- setReadable(GOclusterplot, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

pdf("./results/DE/st/se_all_clusters_dotplot_rna.pdf")
print(dotplot(GOclusterplot))
dev.off()

pdf("./results/DE/st/se_all_clusters_cnet_rna.pdf")
cnetplot(ck)
dev.off()
####se object DE FIN

####ENRICHMENT ANALYSIS DORMANT



