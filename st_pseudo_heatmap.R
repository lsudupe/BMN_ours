




top <- names(sort(selected_genes, decreasing = TRUE))#[1:100]
final_mat_se <- final_mat[top, ] #pseudobulk ordered most var genes
dim(final_mat_se)

## ana heatmap 
matrix <- final_mat[top, ]
ordered_cols <- c("Pg1", "Pg2", "Pg3", "Pg4", "Pg5", "Pg6", "Pg7", "Pg8", "Pg9", "Pg10", "Pg11", "Pg12")
matrix_ordered <- matrix[, ordered_cols]

pheatmap(matrix_ordered, scale = "row", fontsize = 4)

final_mat_se
matrix_ordered

library(factoextra)
library(cluster)
library(ggplot2)

## clusters
a <- fviz_nbclust(matrix_ordered, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(subtitle = "Método del Codo")
b <- fviz_nbclust(matrix_ordered, kmeans, method = "silhouette") +
  labs(subtitle = "Método de la Silueta")
grid.arrange(a, b,ncol = 2)


# Realizar clustering jerárquico
dist_matrix <- dist(t(matrix_ordered))  # Calcula la distancia entre columnas
hc <- hclust(dist_matrix, method = "complete")  # Clustering jerárquico

# Cortar el árbol de clustering para obtener 4 grupos
clusters <- cutree(hc, k = 4)

# Crear el mapa de calor con etiquetas de columnas más grandes y en la parte inferior
a <- pheatmap(matrix_ordered, scale = "row", 
         fontsize_col = 12,  # Ajusta este valor para cambiar el tamaño de las etiquetas de las columnas
         fontsize_row = 4,   # Tamaño de las etiquetas de las filas
         cluster_cols = hc, 
         annotation_col = data.frame(Cluster = as.factor(clusters)),
         show_colnames = TRUE,  # Mostrar nombres de columnas
         show_rownames = TRUE)  # Mostrar nombres de filas

pdf("./results/pseudo/heatmap.pdf", width = 8, height = 15)
print(a)
dev.off()

# Crear el mapa de calor con etiquetas de columnas más grandes y en la parte inferior
pheatmap(matrix_ordered, scale = "row",cutree_cols = 4
            )

