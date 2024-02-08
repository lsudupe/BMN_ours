
## SCRIPT: CytAssist vs Manual

## 08.02.24 Laura Sudupe , git @lsudupe

## Libraries
library(Seurat)
library(tidyverse)
library(STutility)
library(dplyr)
library(ggplot2)


## Define proyect info
ANALYSIS_ID <- "visium_comparation"

## Define project paths
DIR_WD <- getwd()
DIR_ROOT <- file.path(getwd())  
DIR_DATA <- file.path(DIR_ROOT, "data/data/")


## Read data
## samples
samples <- list.files(path = DIR_DATA, recursive = TRUE, pattern = "filtered",full.names = TRUE)
imgs <- list.files(path = DIR_DATA, pattern = "tissue_hires_image.png", recursive = TRUE, full.names = TRUE)
spotfiles <- list.files(path = DIR_DATA, pattern = "tissue_positions_list.csv", recursive = TRUE, full.names = TRUE)
json <- list.files(path = DIR_DATA, pattern = "scalefactors_json.json", recursive = TRUE, full.names = TRUE)
name <- dir(path = DIR_DATA)

## Create the infotable
infoTable <- data.frame(samples = samples, imgs = imgs, spotfiles = spotfiles, json = json, 
                        name = name)

infoTable <- infoTable[is.element(infoTable$name, c("M2_F_2B", "Femur_M2")), ]

## Create the object
se <- InputFromTable(infoTable,
                     platform =  "Visium")

## Add image info
st.object <- GetStaffli(se)
st.object
se <- LoadImages(se, time.resolve = FALSE)

meta <- se@meta.data
meta$name <- gsub("Femur_M2", "CytAssist FFPE", meta$name)
meta$name <- gsub("M2_F_2B", "Manual FFPE", meta$name)

st <- AddMetaData(se, meta)


#QC
ST.FeaturePlot(st, features = "nFeature_RNA", cols = c("lightgray", "mistyrose", "red", "dark red", "black"), ncol = 2, pt.size = 1.3)


p1 <- ggplot() +
  geom_histogram(data = st[[]], aes(nFeature_RNA, color=name), fill = "white", alpha = 0.7, bins = 50) +
  ggtitle("Unique genes per spot")

p2 <- ggplot() +
  geom_histogram(data = st[[]], aes(nCount_RNA, color=name), fill = "white", alpha = 0.7, bins = 50) +
  ggtitle("Total counts per spots")

pdf("./emma_cytvsma/histo_qc.pdf", width = 12, height = 8)
(p1 - p2)
dev.off()

se <- st

# Keep spots with more than 500 unique genes and less than 30% mitochondrial content
se.subset <- SubsetSTData(se, expression = nFeature_RNA > 500 )

cat("Spots removed: ", ncol(se) - ncol(se.subset), "\n")


se <- SCTransform(se)
se <- NormalizeData(se)
se <- ScaleData(se)

#####

manual_data <- subset(se, subset = name == "Manual FFPE")
cytoassist_data <- subset(se, subset = name == "CytAssist FFPE")

# Calculate average expression for each gene in both subsets. Replace `SCT` with the appropriate assay if different.
avg_exp_manual <- AverageExpression(manual_data, return.seurat = FALSE)$SCT
avg_exp_cytoassist <- AverageExpression(cytoassist_data, return.seurat = FALSE)$SCT

# Merge the two datasets to ensure only common genes are analyzed
common_genes <- merge(avg_exp_manual, avg_exp_cytoassist, by = "row.names", all = FALSE)
colnames(common_genes) <- c("gene", "avg_manual", "avg_cytoassist")

# Calculate correlation
correlation <- cor(common_genes$avg_manual, common_genes$avg_cytoassist, method = "pearson")

# Scatter plot of gene expression correlation COMMON GENES
g <- ggplot(common_genes, aes(x = avg_manual, y = avg_cytoassist)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() + scale_y_log10() + # Log scale might be more informative
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  theme_minimal() +
  labs(title = paste("Correlation of Gene Expression: Manual vs. CytoAssist", "\nCorrelation Coefficient:", round(correlation, 3)),
       x = "Average Expression (Manual)",
       y = "Average Expression (CytoAssist)",
       caption = "Each point represents a gene")

pdf("./emma_cytvsma/correlation.pdf", width = 12, height = 8)
print(g)
dev.off()


# Obtener los datos de expresión génica para cada subconjunto
manual_expression <- GetAssayData(manual_data, slot = "counts")
cytoassist_expression <- GetAssayData(cytoassist_data, slot = "counts")

# Contar el número de genes detectados en cada método
# Definir un "gen detectado" como aquel que tiene una expresión mayor que cero en al menos una célula
num_genes_manual <- sum(rowSums(manual_expression > 0) > 0)
num_genes_cytoassist <- sum(rowSums(cytoassist_expression > 0) > 0)


# Imprimir el número de genes
cat("Número de genes detectados en el método Manual:", num_genes_manual, "\n")
cat("Número de genes detectados en el método CytoAssist:", num_genes_cytoassist, "\n")

##comun
# Identificar los genes detectados (expresión mayor que cero) en cada subconjunto
genes_detected_manual <- rownames(manual_expression)[rowSums(manual_expression > 0) > 0]
genes_detected_cytoassist <- rownames(cytoassist_expression)[rowSums(cytoassist_expression > 0) > 0]

# Identificar los genes comunes a ambos métodos
common_genes <- intersect(genes_detected_manual, genes_detected_cytoassist)

# Contar el número de genes comunes
num_common_genes <- length(common_genes)

# Imprimir el número de genes comunes
cat("Número de genes detectados en común entre los métodos Manual y CytoAssist:", num_common_genes, "\n")

    


