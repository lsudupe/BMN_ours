## SCRIPT: Separate hot-spot in graphs BONE MARROW

## 03.12.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(igraph)
library(deldir)
library(dbscan)
library(cluster)
library(colorBlindness)


## Data
mouse <- readRDS("./objects/sp/st_q3.rds")

# Define the subsets to processq
subsets <- c("M1_fem_1C", "M2_F_2B", "M8_F2_1C", "M9_F2_1C")


# Sample
a <- SubsetSTData(mouse, idents = c("M8_F2_1C"))
Idents(object = a) <- a@meta.data[["new_groups"]]
a <- SubsetSTData(a, idents = c("Hot_spot"))

# Extract coordinates
meta <- a@meta.data
coor <- a@tools[["Staffli"]]@meta.data
meta["x_coord"] <- as.vector(coor$x) 
meta["y_coord"] <- as.vector(coor$y) 
coordinates <- as.matrix(meta[, c("x_coord", "y_coord")])

dist_matrix <- as.matrix(dist(coordinates))

# Inicializa una matriz de adyacencia con ceros
adj_matrix <- matrix(0, nrow = nrow(dist_matrix), ncol = nrow(dist_matrix))

# Llena la matriz de adyacencia conectando cada punto solo a sus 6 vecinos más cercanos
for (i in 1:nrow(dist_matrix)) {
  # Ordena y selecciona los índices de los 6 vecinos más cercanos
  neighbors <- order(dist_matrix[i,])[2:7]  # Excluye el propio punto [1]
  adj_matrix[i, neighbors] <- 1
}

# Asegúrate de que la matriz de adyacencia es simétrica
adj_matrix <- pmax(adj_matrix, t(adj_matrix))

# Crea el grafo a partir de la matriz de adyacencia
graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)

# Opcional: Eliminar nodos aislados si los hay
graph <- delete.vertices(graph, which(degree(graph) == 0))

# Visualiza el grafo
plot(graph, vertex.size=5, vertex.label=NA, asp=1)



# Check resolutions
resolutions = seq(0.1, 1, by=0.1)
silhouette_scores = numeric(length(resolutions))

for (j in 1:length(resolutions)) {
  communities <- cluster_louvain(graph, resolution = resolutions[j])
  meta$community <- factor(communities$membership)
  
  # Convertir meta$community en un vector numérico
  community_numeric <- as.numeric(meta$community)
  
  # Calcula la matriz de distancia
  dist_matrix <- dist(coordinates)
  
  # Calcula los coeficientes de silueta
  silhouette_values <- silhouette(community_numeric, dist_matrix)
  
  # Comprobar si hay suficientes clusters para calcular la silueta
  if (length(unique(community_numeric)) > 1) {
    silhouette_scores[j] <- mean(silhouette_values[, "sil_width"])
  } else {
    silhouette_scores[j] <- NA  # NA para resoluciones con un único cluster
  }
}

# Encontrar la resolución con el coeficiente de silueta más alto (ignorando NA)
best_resolution_index <- which.max(silhouette_scores)
best_resolution <- resolutions[best_resolution_index]

# Imprimir o graficar los resultados
plot(resolutions, silhouette_scores, type='b', xlab='Resolution', ylab='Average Silhouette Width')
abline(v = best_resolution, col = "red", lwd = 2)



# Detectar comunidades con el algoritmo de Louvain
communities <- cluster_louvain(graph, resolution = 0.1)
# Visualizar el grafo con las comunidades
plot(communities, graph, vertex.size=5, vertex.label=NA)


# Extract info and add it to your spatial data
community_membership <- communities$membership
community_factor <- as.factor(community_membership)
meta$community <- community_factor
a <- AddMetaData(a, meta)



# Plot spatial
FeatureOverlay(a, features = c("community"), 
               sampleids = 1:6, ncols = 2, 
               pt.size = 1.3)

M9 <- a
# Change level names
# k=15, res=0.2
M8
levels(M8@meta.data$community) <- c("S8_M1", "S8_M2")
# k=15, res=0.5
M2
levels(M2@meta.data$community) <- c("S2_M1", "S2_M2", "S2_M3", "S2_M4", "S2_M5")
# k=15, res=0.2
M1
levels(M1@meta.data$community) <- c("S1_M1", "S1_M2")
# k=15, res=0.5
M9
levels(M9@meta.data$community) <- c("S9_M1", "S9_M2", "S9_M3")

# Save data
saveRDS(M1,"./objects/sp/st_s1_module.rds")
saveRDS(M2,"./objects/sp/st_s2_module.rds")
saveRDS(M8,"./objects/sp/st_s8_module.rds")
saveRDS(M9, "./objects/sp/st_s9_module.rds")



###colorblindness plots
M1 <- graph
communities_M1 <- communities
# Visualizar el grafo con las comunidades
plot(communities_M1, M1, vertex.size=5, vertex.label=NA)


color_mapping <- c(Pg1="#000000", Pg2="#004949", Pg3="#009292", Pg4="#ff6db6", Pg5="#ffb6db",
                   Pg6="#490092", Pg7="#006ddb", Pg8="#b66dff", Pg9="#6db6ff", Pg10="#b6dbff",
                   Pg11="#920000", Pg12="#924900", Pg13="#db6d00", Pg14="#24ff24", Pg15="#ffff6d")


#M1 <- graph
#communities_M1 <- communities
pdf("./results/graphs/communities_M1.pdf")
cols <- color_transparente[membership(communities_M1)]
plot(communities_M1, M1, col = cols,vertex.size=5, mark.col= c("#000000", "#004949"), vertex.label=NA, asp=1, 
     vertex.color = c("#000000", "#006ddb"))
dev.off()

#M2 <- graph
#communities_M2 <- communities
pdf("./results/graphs/communities_M2.pdf")
cols <- color_transparente[membership(communities_M2)]
plot(communities_M2, M2, col = cols,vertex.size=5, mark.col= c("#009292", "#ff6db6", "#ffb6db",
                                                    "#490092", "#006ddb"), vertex.label=NA, asp=1,
     vertex.color = c("#009292", "#ff6db6", "#ffb6db",
                      "#490092", "#006ddb"))
dev.off()

#M8 <- graph
#communities_M8 <- communities
pdf("./results/graphs/communities_M8.pdf")
cols <- color_transparente[membership(communities_M8)]
plot(communities_M8, M8, col = cols,vertex.size=5, mark.col= c("#b66dff","#6db6ff"), vertex.label=NA, asp=1,
     vertex.color = c("#b66dff","#6db6ff"))
dev.off()

#M9 <- graph
#communities_M9 <- communities
pdf("./results/graphs/communities_M9.pdf")
cols <- color_transparente[membership(communities_M9)]
plot(communities_M9, M9, col = cols,vertex.size=5, mark.col= c("#b6dbff","#920000","#924900"), vertex.label=NA, asp=1,
     vertex.color = c("#b6dbff","#920000","#924900"))
dev.off()

#opacidad

color_transparente <- adjustcolor(c("#000000", "#004949"), alpha.f = 0.5)
pdf("./results/graphs/communities_M1_opaco.pdf")
cols <- color_transparente[membership(communities_M1)]
plot(communities_M1, M1, col = cols,vertex.size=5, mark.col= color_transparente, vertex.label=NA, asp=1, 
     vertex.color = NA)
dev.off()

color_transparente <- adjustcolor(c("#009292", "#ff6db6", "#ffb6db",
                                    "#490092", "#006ddb"), alpha.f = 0.5)
pdf("./results/graphs/communities_M2_opaco.pdf")
cols <- color_transparente[membership(communities_M2)]
plot(communities_M2, M2, col = cols,vertex.size=5, mark.col= color_transparente, vertex.label=NA, asp=1,
     vertex.color = NA)
dev.off()

color_transparente <- adjustcolor(c("#b66dff","#6db6ff"), alpha.f = 0.5)
cols <- color_transparente[membership(communities_M8)]
pdf("./results/graphs/communities_M8_opaco.pdf")
plot(communities_M8, M8, col = cols,vertex.size=5, mark.col= color_transparente, vertex.label=NA, asp=1,
     vertex.color = NA)
dev.off()

color_transparente <- adjustcolor(c("#b6dbff","#920000","#924900"), alpha.f = 0.5)
cols <- color_transparente[membership(communities_M9)]
pdf("./results/graphs/communities_M9_opaco.pdf")
plot(communities_M9, M9, col = cols,vertex.size=5, mark.col= color_transparente, vertex.label=NA, asp=1,
     vertex.color = NA)
dev.off()




