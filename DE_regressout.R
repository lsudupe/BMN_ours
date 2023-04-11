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



#Data--------------------------------
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

b <- SetIdent(se, value = se@meta.data[["clustering"]])

###subset the data 
meta_b <- b@meta.data
types_b <- meta_b[,12:21]
types_b["clustering"] <- as.vector(b@meta.data[["clustering"]])

matrix <- as.data.frame(se@assays[["RNA"]]@counts)
matrix_t <- t(matrix)
matrix_t <- as.data.frame(matrix_t)
matrix_t["plasma_value"] <- types_b$MM_MIC

matrix_mini <- as.data.frame(matrix_t[,1:3])
matrix_mini["plasma_value"] <- types_b$MM_MIC

write.csv(matrix_t, "./matrix.csv", row.names=TRUE)

df <- matrix_t
# for value in row and column
# regress out that row final value and 
# new value obtained add it to exact same place in data.frame


# Get the index of the last column (special value)
special_value_col <- ncol(df)

# Define a custom function for the regression
regress_out_special_value <- function(x, i, j, special_value) {
  temp_df <- data.frame(original_value = x, special_value = special_value)
  lm_model <- lm(original_value ~ special_value, data = temp_df)
  
  # Print progress update
  cat(sprintf("Processing row %d, column %d\n", i, j))
  
  return(lm_model$residuals[1])
}

# Iterate through each row and column (except the last one) using the apply function
df[, 1:(special_value_col - 1)] <- apply(df[, 1:(special_value_col - 1)], c(1, 2), function(x, i, j) {
  regress_out_special_value(x, i, j, df[i, special_value_col])
}, i = row(df), j = col(df))

# Print the updated data frame
print(df)












