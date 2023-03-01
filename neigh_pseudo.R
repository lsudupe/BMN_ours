## SCRIPT: Neighbours pseudobulk and DE

## 01.03.23 Laura Sudupe , git @lsudupe

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
se <- readRDS("./objects/sc/integrated/se_deco.rds")


# subset data

se@meta.data[["nbs_2"]] <- factor(se@meta.data[["nbs_2"]], levels = c("NA", "nbs_2", "2"))
vec <- se@meta.data[["nbs_2"]]
levels(vec) <- list(other = "NA", nbs_2 = "nbs_2", cluster = "2")
se@meta.data[["nbs_2"]] <- vec

Seurat::Idents(object = se) <- se@meta.data[["nbs_2"]]
se_sub <- subset(x =se, idents = c("nbs_2", "cluster"))


se_sub@meta.data[["nbs_2"]] <- se_sub@active.ident

# agregate the expression in both
se_sub <- NormalizeData(se_sub)
cct <- AggregateExpression(se_sub, 
                         group.by = c("nbs_2", "name"),
                         assays = 'SCT',
                         slot = "data",
                         return.seurat = FALSE)

# separate matrixes
cts <- cct$SCT
df <- as.data.frame(cts)

# Named vector of sample names
sample <- purrr::set_names(levels(as.factor(se_sub$name)))
sample = as.factor(c(rep(sample,2)))
area <- as.factor(c("nbs_2","nbs_2","nbs_2","nbs_2","cluster","cluster","cluster","cluster"))
mycols = data.frame(row.names = colnames(df),sample, area)
mycols
# Our area of interest
mycols$area <- relevel(mycols$area, "nbs_2")

# Create dsd object
dsd = DESeqDataSetFromMatrix(countData = round(cts),colData = mycols, design = ~ area) 

# Run DESeq2 differential expression analysis
dds <- DESeq(dsd)
#dds$area <- relevel(dds$area, "RZ")

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
colnames(designGO)<-c("nbs_2","cluster")
y <- estimateDisp(y,designGO)

#######RZ######
# contrast of interest
con <- makeContrasts(nbs_2 - cluster,levels=designGO)
fit <- glmQLFit(y,designGO)
qlf <- glmQLFTest(fit,contrast=con)


# nbs_2 - cluster
qlf.nbs_2vscluster <- glmQLFTest(fit,contrast=con[,"nbs_2 - cluster"])
qlf.nbs_2vscluster_matrix <- as.data.frame(topTags(qlf.nbs_2vscluster,n=11330,sort.by = "PValue"))
genes_of_interest_qlf.nbs_2vscluster <- row.names(qlf.nbs_2vscluster_matrix[qlf.nbs_2vscluster_matrix[,"FDR"]<0.005 & qlf.nbs_2vscluster_matrix$logFC>1,])

top.genes  <- unique(genes_of_interest_qlf.nbs_2vscluster)
normalize_matrix <- as.matrix(normalize)
## prepare data
matrix <- normalize_matrix[top.genes, ]

matrix_genes <- data.frame(normalize_matrix) %>%
  rownames_to_column(var = "genes") %>%
  dplyr::filter(genes %in% top.genes)
matrix_genes_ <- data.frame(matrix_genes[,-1], row.names = matrix_genes[,1])

gathered <- matrix_genes %>%
  gather(colnames(matrix_genes)[2:length(colnames(matrix_genes))], key = "samplename", value = "normalize")

##
gathered=within(gathered,{
  area=NA
  area[samplename==c("nbs_2_M1_fem_1C")] = "nbs_2"
  area[samplename==c("nbs_2_M2_F_2B")] = "nbs_2"
  area[samplename==c("nbs_2_M8_F2_1C")] = "nbs_2"
  area[samplename==c("nbs_2_M9_F2_1C")] = "nbs_2"
  area[samplename==c("cluster_M1_fem_1C")] = "cluster"
  area[samplename==c("cluster_M2_F_2B")] = "cluster"
  area[samplename==c("cluster_M8_F2_1C")] = "cluster"
  area[samplename==c("cluster_M9_F2_1C")] = "cluster"
})

# Run pheatmap using the metadata data frame for the annotation
pdf("./results/ST/pseudo/heatmap_neigh2_pseudo.pdf")
print(pheatmap(matrix_genes_, 
               #color = heat_colors, 
               cluster_rows = T, 
               cluster_cols = F,
               show_rownames = T,
               cutree_cols = 4,
               annotation = mycols[, c("area", "sample")], 
               border_color = c("blue"),
               scale = "row", 
               fontsize_row = 5, 
               height = 20)) 
dev.off()


se_sub_st <- SubsetSTData(se, idents = c("nbs_2", "cluster"))

pdf("./results/ST/pseudo/neigh2_genes_M1.pdf")
FeatureOverlay(se_sub_st, features = c("Dram1", "Vill", "Cald1", "Dennd4b"), pt.size = 0.7,  
               sampleids = 1:1,
               ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
dev.off()

pdf("./results/ST/pseudo/neigh2_genes_M2.pdf")
FeatureOverlay(se_sub_st, features = c("Dram1", "Vill", "Cald1", "Dennd4b"), pt.size = 0.7,  
               sampleids = 2:2,
               ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
dev.off()

pdf("./results/ST/pseudo/neigh2_genes_M8.pdf")
FeatureOverlay(se_sub_st, features = c("Dram1", "Vill", "Cald1", "Dennd4b"), pt.size = 0.7,  
               sampleids = 3:3,
               ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
dev.off()

pdf("./results/ST/pseudo/neigh2_genes_M9.pdf")
FeatureOverlay(se_sub_st, features = c("Dram1", "Vill", "Cald1", "Dennd4b"), pt.size = 0.7,  
               sampleids = 4:4,
               ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
dev.off()


#####ALL samples
se <- readRDS("./objects/sc/integrated/se_deco.rds")

pdf("./results/ST/pseudo/neigh2_genes_M1_all.pdf")
FeatureOverlay(se, features = c("Dram1", "Vill", "Cald1", "Dennd4b"), pt.size = 0.7,  
               sampleids = 1:1,
               ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
dev.off()

pdf("./results/ST/pseudo/neigh2_genes_M2_all.pdf")
FeatureOverlay(se, features = c("Dram1", "Vill", "Cald1", "Dennd4b"), pt.size = 0.7,  
               sampleids = 2:2,
               ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
dev.off()

pdf("./results/ST/pseudo/neigh2_genes_M3_F_all.pdf")
FeatureOverlay(se, features = c("Dram1", "Vill", "Cald1", "Dennd4b"), pt.size = 0.7,  
               sampleids = 3:3,
               ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
dev.off()

pdf("./results/ST/pseudo/neigh2_genes_M3_fem_all.pdf")
FeatureOverlay(se, features = c("Dram1", "Vill", "Cald1", "Dennd4b"), pt.size = 0.7,  
               sampleids = 4:4,
               ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
dev.off()

pdf("./results/ST/pseudo/neigh2_genes_M8_all.pdf")
FeatureOverlay(se, features = c("Dram1", "Vill", "Cald1", "Dennd4b"), pt.size = 0.7,  
               sampleids = 5:5,
               ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
dev.off()

pdf("./results/ST/pseudo/neigh2_genes_M9_all.pdf")
FeatureOverlay(se, features = c("Dram1", "Vill", "Cald1", "Dennd4b"), pt.size = 0.7,  
               sampleids = 6:6,
               ncols = 3, cols = c("darkblue", "cyan", "yellow", "red", "darkred"))
dev.off()

