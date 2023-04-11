## SCRIPT: Plot gradients results BM project

## 19.03.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(dplyr)
library(hrbrthemes)
library(extrafont)
library(ggthemes)
library(viridis)


#Data---------------------------------
se <- readRDS("./objects/sc/integrated/se_deco.rds")
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

########cell types driving the cluster
x <- se

##replace all values below 10% with 0
meta <- x@meta.data
a <- meta[,12:21]
a["cluster"] <- as.vector(x@meta.data[["clustering"]])

df <- a
# Filter values below 10%

# Calculate the 75th percentile threshold for each row
thresholds <- apply(df[, -ncol(df)], 1, function(x) quantile(x, 0.75))

# Filter values below the 75th percentile threshold for each row
df_filtered <- df %>%
  mutate(across(-cluster, ~ifelse(. < thresholds[row_number()], 0, .)))

# Filter values below 10%
df_filtered <- df %>%
  mutate(across(-cluster, ~ifelse(. < 0.15, 0, .)))

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
pdf(file.path("./results/ST/gradient/cell_types_per_cluster_15percent.pdf"))
print(plot)
dev.off()


########


library(RColorBrewer)
display.brewer.all()
color <- brewer.pal(11,"Spectral")
color <- rev(color)

p <-FeatureOverlay(se, features = c("MM_MIC"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = color)
p <-FeatureOverlay(se, features = c("MM_MIC"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = c("lightgray", "mistyrose", "red", "dark red", "black"))

p <- ST.FeaturePlot(se, features = c("MM_MIC"),indices  = 1:6, pt.size = 0.80,ncol = 2 , grid.ncol = 1, palette = "Spectral")
p <- ST.FeaturePlot(se, features = c("MM_MIC"),indices  = 1:6, pt.size = 0.80,ncol = 2 , grid.ncol = 1, cols = c("lightgray", "mistyrose", "red", "dark red", "black"))



pdf(file.path("./results/ST/gradient/gradient_rojo.pdf"))
print(p)
dev.off()

####marker genes https://www.nature.com/articles/s41467-022-33944-z
#CD81, XBP1, MKI67, PCNA, FLNA, CD44
###ANORMAL plasma cells
p <-FeatureOverlay(se, features = c("Xbp1"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = color)

pdf(file.path("./results/ST/gradient/pc_markers/Xbp1.pdf"))
print(p)
dev.off()

p <-FeatureOverlay(se, features = c("Cd81"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = color)

pdf(file.path("./results/ST/gradient/pc_markers/Cd81.pdf"))
print(p)
dev.off()

p <-FeatureOverlay(se, features = c("Mki67"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = color)

pdf(file.path("./results/ST/gradient/pc_markers/Mki67.pdf"))
print(p)
dev.off()

p <-FeatureOverlay(se, features = c("Pcna"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = color)

pdf(file.path("./results/ST/gradient/pc_markers/Pcna.pdf"))
print(p)
dev.off()

p <-FeatureOverlay(se, features = c("Flna"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = color)

pdf(file.path("./results/ST/gradient/pc_markers/Flna.pdf"))
print(p)
dev.off()

p <-FeatureOverlay(se, features = c("Cd44"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = color)

pdf(file.path("./results/ST/gradient/pc_markers/Cd44.pdf"))
print(p)
dev.off()


####marker genes FIN

meta <- se@meta.data
meta_mm <- meta[grepl("MM", meta$condition),]
meta_healthy <- meta[grepl("control", meta$condition),]

## density plots
# Make the histogram
M1_fem_1C <- meta[grepl("M1_fem_1C", meta$name),]
M2_F_2B <- meta[grepl("M2_F_2B", meta$name),]
M3_F_1C <- meta[grepl("M3_F_1C", meta$name),]
M3_fem_1C <- meta[grepl("M3_fem_1C", meta$name),]
M8_F2_1C <- meta[grepl("M8_F2_1C", meta$name),]
M9_F2_1C <- meta[grepl("M9_F2_1C", meta$name),]

pdf(file.path("./results/ST/gradient/density_together.pdf"))
print(ggplot(data=meta, aes(x=MM_MIC, group=name, fill=name)) +
  geom_density(adjust=1.5, alpha=.4) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()+
  ylab("") +
  xlab("MM values (%)"))
dev.off()

# Using Small multiple
pdf(file.path("./results/ST/gradient/density_mm.pdf"))
ggplot(data=meta_mm, aes(x=MM_MIC, group=name, fill=name)) +
  geom_density(adjust=1.5) +
  theme_ipsum() +
  facet_wrap(~name) +
  xlim(0, 1) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank()
  )
dev.off()

##add list values to df
for (i in 1:length(lista)){
  a <- lista[[i]]
  pdf(paste("./results/ST/gradient/", names(lista[i]),"_density.pdf",sep=""))
  a %>%
    ggplot( aes(x=MM_MIC)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
    ylab("") +
    xlab("MM values (%)")
  dev.off()
}


ex

