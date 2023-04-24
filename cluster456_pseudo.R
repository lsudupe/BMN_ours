## SCRIPT: Cluster 4vs5vs6 pseudobulk Bone Marrow project

## 23.04.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(edgeR)
library(dplyr)
library(data.table)
library(sva)
library(STutility)

#Data---------------------------------
all <- readRDS("./objects/heterogeneity/se_hierarchical.rds")
b <- SetIdent(all, value = all@meta.data[["clustering"]])
b <- SubsetSTData(b, idents = c("4", "5","6"))
b <- SetIdent(b, value = b@meta.data[["name"]])
subset <- SubsetSTData(b, idents = c("M2_F_2B","M8_F2_1C"))
DefaultAssay(subset) <- "RNA"
#subset <- ScaleData(subset)
subset <- SetIdent(subset, value = subset@meta.data[["clustering"]])

subset@meta.data[["pseudo"]] <- subset@active.ident
vec <- subset@meta.data[["pseudo"]]
levels(vec) <- list(cluster4 = "4", cluster5 = "5", cluster6 = "6")
subset@meta.data[["pseudo"]] <- vec
Idents(subset) <-subset@meta.data[["pseudo"]]

# agregate the expression 
cct_s <- AggregateExpression(subset, 
                           assays = 'RNA',
                           group.by = c("name", "pseudo"),
                           slot = "count",
                           return.seurat = FALSE)
# separate matrixes
cts <- cct_s$RNA
df <- as.data.frame(cts)

# Named vector of sample names
sample <- purrr::set_names(levels(as.factor(subset$name)))
sample = as.factor(c(rep(sample,3)))
area <- as.factor(c("cluster4","cluster4","cluster5","cluster5","cluster6","cluster6"))
mycols = data.frame(row.names = colnames(df),sample, area)
mycols

# correct columns
mycols$sample <- c("M2_F_2B","M2_F_2B","M2_F_2B","M8_F2_1C","M8_F2_1C","M8_F2_1C")
mycols$area <- c("cluster4","cluster5","cluster6","cluster4","cluster5","cluster6")
mycols

# Create dsd object
dsd = DESeqDataSetFromMatrix(countData = cts, colData = mycols, design = ~ area) 

# Run DESeq2 differential expression analysis
dds <- DESeq(dsd)
dataset <- counts(dds, 
                  normalized = TRUE)

## filter genes
cpmgo <- cpm(dataset)
numbermorezero <- apply(cpmgo,1,function(x){sum(x>1)})
hist(numbermorezero,breaks=200)  
dataset <- dataset[numbermorezero>2,]  
hist(dataset,breaks=200)  

## batch correction
adjusted_counts <- ComBat_seq(dataset, batch=mycols$sample,
                              group=mycols$area)
# Normalize the counts
normalize <- NormalizeData(adjusted_counts)

## EDGER
y <- DGEList(counts=adjusted_counts,group=mycols$area)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
designGO <- model.matrix(~ 0+ mycols$area)
colnames(designGO)<-c("cluster4","cluster5","cluster6")
y <- estimateDisp(y,designGO)

#######cluster4######
# contrast of interest
con <- makeContrasts(cluster4 - cluster5, cluster4 - cluster6, cluster5 - cluster6,levels=designGO)
fit <- glmQLFit(y,designGO)
qlf <- glmQLFTest(fit,contrast=con)
# 4-5
qlf.4vs5 <- glmQLFTest(fit,contrast=con[,"cluster4 - cluster5"])
qlf.4vs5_matrix <- as.data.frame(topTags(qlf.4vs5,n=11336,sort.by = "PValue"))
genes_of_interest_4vs5 <- row.names(qlf.4vs5_matrix[qlf.4vs5_matrix[,"PValue"]<0.0005 & qlf.4vs5_matrix$logFC>1,])
# 4-6
qlf.4vs6 <- glmQLFTest(fit,contrast=con[,"cluster4 - cluster6"])
qlf.4vs6_matrix <- as.data.frame(topTags(qlf.4vs6,n=11336,sort.by = "PValue"))
genes_of_interest_4vs6 <- row.names(qlf.4vs6_matrix[qlf.4vs6_matrix[,"PValue"]<0.0005 & qlf.4vs6_matrix$logFC>1,])
#
top.genes  <- c(genes_of_interest_4vs5, genes_of_interest_4vs6)
top.genes  <- unique(top.genes)
#######cluster4###### FIN

#######cluster5######
# contrast of interest
con <- makeContrasts(cluster5 - cluster4, cluster5 - cluster6, cluster4 - cluster6,levels=designGO)
fit <- glmQLFit(y,designGO)
qlf <- glmQLFTest(fit,contrast=con)
# 5-4
qlf.5vs4 <- glmQLFTest(fit,contrast=con[,"cluster5 - cluster4"])
qlf.5vs4_matrix <- as.data.frame(topTags(qlf.5vs4,n=11336,sort.by = "PValue"))
genes_of_interest_5vs4 <- row.names(qlf.5vs4_matrix[qlf.5vs4_matrix[,"PValue"]<0.0005 & qlf.5vs4_matrix$logFC>1,])
# 5-6
qlf.5vs6 <- glmQLFTest(fit,contrast=con[,"cluster5 - cluster6"])
qlf.5vs6_matrix <- as.data.frame(topTags(qlf.5vs6,n=11336,sort.by = "PValue"))
genes_of_interest_5vs6 <- row.names(qlf.5vs6_matrix[qlf.5vs6_matrix[,"PValue"]<0.0005 & qlf.5vs6_matrix$logFC>1,])
#
top.genes  <- c(genes_of_interest_5vs4, genes_of_interest_5vs6)
top.genes  <- unique(top.genes)
#######cluster5###### FIN

#######cluster6######
# contrast of interest
con <- makeContrasts(cluster6 - cluster5, cluster6 - cluster4, cluster5 - cluster4,levels=designGO)
fit <- glmQLFit(y,designGO)
qlf <- glmQLFTest(fit,contrast=con)
# 6-5
qlf.6vs5 <- glmQLFTest(fit,contrast=con[,"cluster6 - cluster5"])
qlf.6vs5_matrix <- as.data.frame(topTags(qlf.6vs5,n=11336,sort.by = "PValue"))
genes_of_interest_6vs5 <- row.names(qlf.6vs5_matrix[qlf.6vs5_matrix[,"PValue"]<0.0005 & qlf.6vs5_matrix$logFC>1,])
# 6-4
qlf.6vs4 <- glmQLFTest(fit,contrast=con[,"cluster6 - cluster4"])
qlf.6vs4_matrix <- as.data.frame(topTags(qlf.6vs4,n=11336,sort.by = "PValue"))
genes_of_interest_6vs4 <- row.names(qlf.6vs4_matrix[qlf.6vs4_matrix[,"PValue"]<0.0005 & qlf.6vs4_matrix$logFC>1,])
#
top.genes  <- c(genes_of_interest_6vs5, genes_of_interest_6vs4)
top.genes  <- unique(top.genes)
#######cluster6###### FIN


## prepare data
normalize_matrix <- as.matrix(normalize)
matrix <- normalize_matrix[top.genes, ]

matrix_genes <- data.frame(normalize_matrix) %>%
  rownames_to_column(var = "genes") %>%
  dplyr::filter(genes %in% top.genes)

## remove 0 values 
matrix_genes <- matrix_genes[matrix_genes$M2_F_2B_cluster4 != 0 & 
                               matrix_genes$M2_F_2B_cluster5 != 0 & 
                               matrix_genes$M8_F2_1C_cluster4 != 0 & 
                               matrix_genes$M8_F2_1C_cluster5 != 0 & 
                               matrix_genes$M2_F_2B_cluster6 != 0 & 
                               matrix_genes$M8_F2_1C_cluster6 != 0, ]

matrix_genes_ <- data.frame(matrix_genes[,-1], row.names = matrix_genes[,1])

my_col_order <- c("M2_F_2B_cluster4","M8_F2_1C_cluster4", "M2_F_2B_cluster5","M8_F2_1C_cluster5",
                  "M2_F_2B_cluster6","M8_F2_1C_cluster6")

# Run pheatmap using the metadata data frame for the annotation
pdf("./prueba_6.pdf")
print(pheatmap(matrix_genes_[my_col_order], 
               #color = heat_colors, 
               cluster_rows = T, 
               cluster_cols = F,
               show_rownames = T,
               cutree_cols = 4,
               annotation = mycols[, c("area", "sample")],
               scale = "row", 
               fontsize_row = 5, 
               height = 20)) 
dev.off()
