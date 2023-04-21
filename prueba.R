


# Create a sample dataframe
set.seed(42)
genes <- 10
spots <- 5
data <- as.data.frame(matrix(runif((genes + 1) * spots), ncol = spots))
colnames(data) <- paste0("Spot", 1:spots)
rownames(data) <- c(paste0("Gene", 1:genes), "Plasma_Percentage")
####EXAMPLE DATA FIN

# Convert the data frame to a matrix
data_matrix <- as.matrix(data)
b <- data_matrix
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

regress <- CreateAssayObject(data)
b@assays[["regress"]] <- regress
b@assays$regress@key <- "regress_"

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

## separate list
list2env(objects,envir=.GlobalEnv)
Idents(object = M2_F_2B)  <- M2_F_2B@meta.data[["clustering"]]         
DefaultAssay(M2_F_2B) <- "regress"
DefaultAssay(M2_F_2B) <- "SCT"

de.subset <- SubsetSTData(M2_F_2B, idents = c("4", "5","6"))

markers <- FindAllMarkers(de.subset,min.pct = 0.1, logfc.threshold = 0.25)

markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
cluster <- subset(markers, p_val_adj < 0.05 & 0.05 < avg_log2FC)
top100 <- cluster %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)


##################################################################
#Subsetting top 100 markers with adjusted p values lower than .05#
##################################################################
top100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)

df <- top100pval[,7:6]
dfsample <- split(df$gene,df$cluster)
length(dfsample)

DefaultAssay(de.subset) <- "SCT"

pdf("./prueba.pdf")
DoHeatmap(de.subset, features = top100pval$gene)#,disp.min = -2, disp.max = 2)
dev.off()

#The output of length(dfsample) returns how many clusters you have
#Here there at 9 clusters (0, 1, 2, 3, 4, 5, 6, 7 and 8)
#I'm sure there's a better way but you have to make a line like below for each cluster
library(clusterProfiler)
library(org.Mm.eg.db)

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
pdf("./prueba_2.pdf")
print(dotplot(GOclusterplot))
dev.off()

pdf("./prueba3.pdf")
cnetplot(ck)
dev.off()

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
dotplot(KEGGclusterplot)

#########clusterprofile FIN

