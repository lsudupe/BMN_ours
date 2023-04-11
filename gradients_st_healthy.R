## SCRIPT: Plot gradients results BM project healthy samples

## 02.04.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(dplyr)
library(hrbrthemes)
library(extrafont)
library(ggthemes)
library(viridis)
library("scales")

#Data---------------------------------
se <- readRDS("./objects/heterogeneity/healthy/se_hierarchical.rds")
x <- se

##replace all values below 10% with 0
meta <- x@meta.data
a <- meta[,12:20]
a["cluster"] <- as.vector(x@meta.data[["clustering"]])

####which cell types are driven my clusters

df <- a
# Filter values below 10%

# Calculate the 75th percentile threshold for each row
thresholds <- apply(df[, -ncol(df)], 1, function(x) quantile(x, 0.95))

# Filter values below the 75th percentile threshold for each row
df_filtered <- df %>%
  mutate(across(-cluster, ~ifelse(. < thresholds[row(.)], 0, .)))

# Reshape the data frame to a long format for plotting
df_long <- df_filtered %>%
  pivot_longer(cols = -cluster,
               names_to = "cell_type",
               values_to = "percentage") %>%
  filter(percentage > 0)

# Create a plot to visualize the number of major cell types driving each cluster
plot <- ggplot(df_long, aes(x = cluster, y = cell_type)) +
  geom_count(aes(size = ..n.., color = cell_type)) +
  labs(title = "Major Cell Types Driving Each Cluster",
       x = "Cluster",
       y = "Cell Type",
       size = "Count") +
  theme_minimal() +
  scale_color_discrete(name = "Cell Type")

# Print the plot
pdf(file.path("./results/ST/gradient/healthy/cell_types_per_cluster_95_quantile.pdf"))
print(plot)
dev.off()

####

a[a < 0.1] <- 0
b <- apply(a, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

meta[,12:20] <- b
x@meta.data <- meta

#x.subset <- SubsetSTData(x, expression = Tcell >= 0.1)
#se.subset <- SubsetSTData(se, expression = nFeature_RNA >= 2000)

####cor
matrix <- data.matrix(b, rownames.force = NA)
M <- cor(matrix)

pdf(file.path("./results/ST/gradient/healthy/prueba_cor.pdf"))
print(corrplot(M, method = 'square', title="Cell type correlation",
               col=colorRampPalette(c("blue","white","red"))(100),
               order = 'original', type = 'lower', diag = FALSE,
               mar=c(0,0,1,0)))
dev.off()

####cor fin

####spatial plots


pdf(file.path("./results/ST/gradient/healthy/prueba.pdf"))
print(ST.FeaturePlot(x, features = c("Erythroblasts", "Neutrophils"),sampleids = 1, pt.size = 0.7,ncol = 1 , 
                     value.scale = "all" , palette = "Spectral"))
dev.off()

p1 <- ST.FeaturePlot(x, features = c("Erythroblasts", "Neutrophils"),sampleids = 1, pt.size = 0.5,ncol = 1 , 
                     value.scale = "all" , palette = "Spectral")
p2 <- VlnPlot(x, features = c("Erythroblasts", "Neutrophils"), ncol = 2, group.by = "clustering")

pdf(file.path("./results/ST/gradient/healthy/violin_plot_celltypes_per_cluster.pdf"))
VlnPlot(x, features = c("Tcell", "Bcell","Erythroblasts" ,"Monocytes","Neutrophils", 
                        "MSC","EC", "DC" ,"NK"), ncol = 3, group.by = "clustering")
dev.off()

library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(11,"Spectral")
color <- rev(color)

pdf(file.path("./results/ST/gradient/healthy/gradiente_neutro_1.pdf"))
print(FeatureOverlay(x, features = c("Neutrophils"), sampleids = 1,pt.size = 1.5,ncol = 1 , 
                     value.scale = "all" ,cols = color))
dev.off()

pdf(file.path("./results/ST/gradient/healthy/gradiente_neutro_2.pdf"))
print(FeatureOverlay(x, features = c("Neutrophils"), sampleids = 2,pt.size = 1.5,ncol = 1 , 
                     value.scale = "all" ,cols = color))
dev.off()

pdf(file.path("./results/ST/gradient/healthy/gradiente_Erythroblasts_2.pdf"))
print(FeatureOverlay(x, features = c("Erythroblasts"), sampleids = 2,pt.size = 1.5,ncol = 1 , 
                     value.scale = "all" ,cols = color))
dev.off()

pdf(file.path("./results/ST/gradient/healthy/gradiente_Bcell_2.pdf"))
print(FeatureOverlay(x, features = c("MSC"), sampleids = 2,pt.size = 1.5,ncol = 1 , 
                     value.scale = "all" ,cols = color))
dev.off()

pdf(file.path("./results/ST/gradient/healthy/gradiente_Monocytes_2.pdf"))
print(FeatureOverlay(x, features = c("Monocytes"), sampleids = 2,pt.size = 1.5,ncol = 1 , 
                     value.scale = "all" ,cols = color))
dev.off()

HSVPlot(x, features = c("Neutrophils","Erythroblasts"), sampleids = 1:2, dark.theme = FALSE, 
        add.alpha = TRUE, show.sb=TRUE)

pdf(file.path("./results/ST/gradient/healthy/prueba_3.pdf"))
print(p1 - p2)
dev.off()
