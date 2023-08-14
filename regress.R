## SCRIPT: Regress out Bone Marrow project

## 14.08.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(STutility)
library(RColorBrewer)
library(UCell)

# Data
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")
x <- se

DefaultAssay(x) <- "SCT"
x <- ScaleData(x)
###subset the data 
meta <- x@meta.data
types <- meta[,12:18]
types["clustering"] <- as.vector(x@meta.data[["clustering"]])


###REGRESS OUT PCs
###extract counts
matrix <- as.data.frame(x@assays[["SCT"]]@counts)
matrix_t <- t(matrix)
matrix_t <- as.data.frame(matrix_t)
matrix_t["plasma_value"] <- types$MM_MIC
matrix_final <- as.data.frame(t(matrix_t))

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
### add new infor to the object
regress <- CreateAssayObject(data)
x@assays[["regress_pc"]] <- regress
x@assays$regress_pc@key <- "regress_pc_"

###REGRESS OUT PCs FIN

###REGRESS OUT monocytes
###extract counts
matrix <- as.data.frame(x@assays[["SCT"]]@counts)
matrix_t <- t(matrix)
matrix_t <- as.data.frame(matrix_t)
matrix_t["monocytes"] <- types$Monocytes
matrix_final <- as.data.frame(t(matrix_t))

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
### add new infor to the object
regress <- CreateAssayObject(data)
x@assays[["regress_monocytes"]] <- regress
x@assays$regress_monocytes@key <- "regress_monocytes_"
###REGRESS OUT monocytes FIN

##extract m8
#M1
Idents(object = x) <- x@meta.data[["name"]]
m8 <- SubsetSTData(object = x, ident = c("M8_F2_1C"))

##plotting
DefaultAssay(m8) <- "SCT"
p <- VlnPlot(m8, features = c("Mmp9"), group.by = "clustering")
pdf(paste("./results/regress/violin_Mmp9_normal.pdf",sep=""))
print(p)
dev.off()

p <- VlnPlot(m8, features = c("Cd44"), group.by = "clustering")
pdf(paste("./results/regress/violin_Cd44_normal.pdf",sep=""))
print(p)
dev.off()

DefaultAssay(m8) <- "regress_pc"
p <- VlnPlot(m8, features = c("Cd44"), group.by = "clustering")
pdf(paste("./results/regress/violin_Cd44_regress.pdf",sep=""))
print(p)
dev.off()

DefaultAssay(m8) <- "regress_monocytes"
p <- VlnPlot(m8, features = c("Mmp9"), group.by = "clustering")
pdf(paste("./results/regress/violin_Mmp9_regress.pdf",sep=""))
print(p)
dev.off()

