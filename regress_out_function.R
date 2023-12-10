## SCRIPT: Regress out variable from counts

## 23.11.23 Laura Sudupe , git @lsudupe

regress_out <- function(a, variable_column_name) {
  # Extract counts
  matrix <- as.data.frame(a@assays[["RNA"]]@counts)
  matrix <- ScaleData(matrix)
  matrix_t <- t(matrix)
  matrix_t <- as.data.frame(matrix_t)
  matrix_t["variable_value"] <- a@meta.data[[variable_column_name]]
  matrix_final <- as.data.frame(t(matrix_t))
  
  # Regress out
  data <- matrix_final
  # Get the number of genes and spots
  genes <- nrow(data) - 1 # Subtract 1 to exclude the last row with variable values
  spots <- ncol(data)
  # Convert the data frame to a matrix
  data_matrix <- as.matrix(data)
  # Initialize the residuals matrix with the same dimensions as the data matrix (without the last row)
  residuals_matrix <- matrix(0, nrow = genes, ncol = spots)
  # Loop through each gene
  for (i in 1:genes) {
    # Create a linear model to regress out the specified variable
    model <- lm(data_matrix[i, 1:spots] ~ data_matrix[genes + 1, 1:spots])
    # Save the residuals in the corresponding row of the residuals_matrix
    residuals_matrix[i, ] <- model$residuals
  }
  # Replace the original data_matrix with the residuals
  data_matrix[1:genes, 1:spots] <- residuals_matrix
  # Convert the matrix back to a data frame
  data <- as.data.frame(data_matrix)
  
  # Add new info to the object
  regress <- CreateAssayObject(data)
  a@assays[["regress"]] <- regress
  a@assays$regress@key <- "regress_"
  # Scale
  DefaultAssay(a) <- "regress"
  
  return(a)
}
