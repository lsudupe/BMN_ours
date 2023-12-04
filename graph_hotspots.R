## SCRIPT: Separate hot-spot in graphs BONE MARROW

## 03.12.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(igraph)
library(deldir)
library(dbscan)


## Data
mouse <- readRDS("./objects/sp/st_q3.rds")

# Define the subsets to processq
subsets <- c("M1_fem_1C", "M2_F_2B", "M8_F2_1C", "M9_F2_1C")

# Sample
M8 <- SubsetSTData(mouse, idents = c("M8_F2_1C"))
Idents(object = M8) <- M8@meta.data[["new_groups"]]
M8 <- SubsetSTData(M8, idents = c("Hot_spot"))

# Extract coordinates
meta <- M8@meta.data
coor <- M8@tools[["Staffli"]]@meta.data
meta["x_coord"] <- as.vector(coor$x) 
meta["y_coord"] <- as.vector(coor$y) 
coordinates <- as.matrix(meta[, c("x_coord", "y_coord")])

# Define the number of nearest neighbors 
k <- 15  

# Find k nearest neighbors for each spot
# The function returns a list of integer vectors, each containing the indices of the k nearest neighbors
knn_result <- kNN(coordinates, k=k, search="kd")

# Initialize an empty graph
graph <- make_empty_graph(nrow(coordinates), directed = FALSE)

# Add vertices with spatial information as attributes
V(graph)$x <- meta$x_coord
V(graph)$y <- meta$y_coord

# Iterate through the knn_result list to add edges
for (i in 1:length(knn_result)) {
  # Get the indices of the k-nearest neighbors for the i-th spot
  # Exclude the first one because it's the point itself
  neighbors <- knn_result[[i]][-1]
  # Add edges from the current point to its neighbors
  graph <- graph + edges(cbind(rep(i, length(neighbors)), neighbors))
}

# Simplify the graph to remove loops and multiple edges
graph <- simplify(graph)

# Perform community detection using the Louvain algorithm
communities <- cluster_louvain(graph,resolution =  0.2)

# Add the community information to the metadata
meta$community <- factor(communities$membership)

# Extract the layout from the graph's vertex attributes
layout_matrix <- cbind(V(graph)$x, V(graph)$y)

# Plot the graph using the layout based on the original spatial coordinates
#plot(graph, layout=layout_matrix, vertex.size=3, vertex.label=NA, edge.arrow.size=.1)
# Color the vertices based on the community membership to visualize clusters
plot(communities, graph, layout=layout_matrix, vertex.size=3, vertex.label=NA, edge.arrow.size=.1, edge.width=0.1)

# Extract info and add it to your spatial data
community_membership <- communities$membership
community_factor <- as.factor(community_membership)
meta$community <- community_factor
M8 <- AddMetaData(M8, meta)

# Plot spatial
FeatureOverlay(M8, features = c("community"), 
               sampleids = 1:6, ncols = 2, 
               pt.size = 1.3)


# Change level names
# k=15, res=0.2
M8
levels(M8@meta.data$community) <- c("S8_M1", "S8_M2")
# k=15, res=0.5
M2
levels(M2@meta.data$community) <- c("S2_M1", "S2_M2", "S2_M3", "S2_M4")
# k=15, res=0.2
M1
levels(M1@meta.data$community)[levels(M1@meta.data$community) == "1"] <- "S1_M1"
# k=15, res=0.5
M9
levels(M9@meta.data$community) <- c("S9_M1", "S9_M2", "S9_M3")

# Save data
saveRDS(M1,"./objects/sp/st_s1_module.rds")
saveRDS(M2,"./objects/sp/st_s2_module.rds")
saveRDS(M8,"./objects/sp/st_s8_module.rds")
saveRDS(M9, "./objects/sp/st_s9_module.rds")

