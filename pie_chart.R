## SCRIPT: Piechart results BM project

## 13.04.23 Laura Sudupe , git @lsudupe

#Libraries---------------------------------
library(Seurat)
library(ggplot2)
library(STutility)
library(tidyverse)

#Data---------------------------------
se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")
x <- se

###coordinates
coor <- se@tools[["Staffli"]]@meta.data

##replace all values below 10% with 0
meta <- x@meta.data
a <- meta[,12:18]
a["cluster"] <- as.vector(x@meta.data[["clustering"]])

df <- a

############ Filter values below whatever
# Calculate the 75th percentile threshold for each row
#thresholds <- apply(df[, -ncol(df)], 1, function(x) quantile(x, 0.75))

# Filter values below the 75th percentile threshold for each row
#df_filtered <- df %>%
 # mutate(across(-cluster, ~ifelse(. < thresholds[row_number()], 0, .)))

df["x_coord"] <- as.vector(coor$x) 
df["y_coord"] <- as.vector(coor$y) 
df["sample"] <- as.vector(meta$name) 

# Filter values below the 15% threshold for each row
df_filtered <- df %>%
  mutate(across(-c(cluster, x_coord, y_coord, sample), ~ifelse(. < 0.20, 0, .)))

# Rescale non-zero values in each row to sum up to 1 and overwrite the original values
df_rescaled <- df_filtered %>%
  rowwise() %>%
  mutate(across(-c(cluster, x_coord, y_coord, sample), 
                ~ifelse(. != 0, . / sum(c_across(-c(cluster, x_coord, y_coord, sample))[c_across(-c(cluster, x_coord, y_coord, sample)) != 0]), 0))) %>%
  ungroup()

# Reshape the data frame to a long format for plotting
df_new <- df_rescaled
df_new$x_coord <- NULL
df_new$y_coord <- NULL
df_new$sample <- NULL

df_long <- df_new %>%
  pivot_longer(cols = -cluster,
               names_to = "cell_type",
               values_to = "percentage") %>%
  filter(percentage > 0)

##color
cells_order <- c("Bcell", "DC", "Erythroblasts","MM_MIC", "Monocytes","Neutrophils","Tcell")
cell_type_colors <- brewer.pal(length(cells_order), "RdBu")
cell_type_colors <-  c("#ef8a62" ,"#ffd1b6" ,"#faeae0" ,"#b2182b" ,"#d1e5f0" ,"#67a9cf" ,"#2166ac")
cell_type_color_map <- setNames(cell_type_colors, cells_order)

# Create a plot to visualize the number of major cell types driving each cluster
plot <- ggplot(df_long, aes(x = cluster, y = cell_type, color = cell_type)) +
  geom_count(aes(size = ..n..)) +
  labs(title = "Major Cell Types Driving Each Cluster",
       x = "Cluster",
       y = "Cell Type",
       size = "Count") +
  theme_minimal() +
  scale_color_manual(name = "Cell Type", values = cell_type_color_map)

# Print the plot
pdf(file.path("./results/ST/gradient/cell_types_per_cluster_20percent_color.pdf"))
print(plot)
dev.off()

# Separate the dataframe by sample
df_list <- split(df_rescaled, df_rescaled$sample)

# Access the dataframe for each sample by its name
df_M1_fem_1C <- df_list[["M1_fem_1C"]]
df_M2_F_2B <- df_list[["M2_F_2B"]]
df_M3_F_1C <- df_list[["M3_F_1C"]]
df_M3_fem_1C <- df_list[["M3_fem_1C"]]
df_M8_F2_1C <- df_list[["M8_F2_1C"]]
df_M9_F2_1C <- df_list[["M9_F2_1C"]]

########################
for (i in 1:length(df_list)){
  data <- df_list[[i]]

  # Reshape the data frame to a long format for plotting
  
  data$sample <- NULL
  data$cluster <- NULL
  colnames(data)[which(names(data) == "x_coord")] <- "x"
  colnames(data)[which(names(data) == "y_coord")] <- "y"

  num_cell_types <- ncol(data) - 2

  data_long <- data %>% 
    tidyr::pivot_longer(cols = 1:num_cell_types, names_to = "cell_type", values_to = "percentage")

  data_polar <- data_long %>%
    group_by(x, y) %>%
    mutate(end = cumsum(2 * pi * percentage),
          start = lag(end, default = 0)) %>%
    ungroup()

  # Filter out values equal to or below 0
  data_polar <- data_polar %>%
    filter(percentage > 0)

  num_points <- 100
  radius <- 0.73 

  data_pie <- tidyr::expand_grid(x = unique(data_polar$x),
                                y = unique(data_polar$y),
                                 point = seq_len(num_points + 1)) %>%
   left_join(data_polar, by = c("x", "y")) %>%
   mutate(angle = ifelse(point == 1, 0, start + (end - start) * (point - 2) / (num_points - 1)),
          x_plot = x + radius * cos(angle) * (point != 1),
           y_plot = y + radius * sin(angle) * (point != 1))

  #library(RColorBrewer)
  #cell_type_colors <- brewer.pal(num_cell_types, "RdBu")
  cells_order <- c("Bcell", "DC", "Erythroblasts","MM_MIC", "Monocytes","Neutrophils","Tcell")
  #cell_type_color_map <- setNames(cell_type_colors, cells_order)
  cell_type_colors <-  c("#ef8a62" ,"#ffd1b6" ,"#faeae0" ,"#b2182b" ,"#d1e5f0" ,"#67a9cf" ,"#2166ac")
  cell_type_color_map <- setNames(cell_type_colors, cells_order)

  pie_chart_plot <- ggplot(data_pie, aes(x = x_plot, y = y_plot, group = interaction(x, y, cell_type))) +
    geom_polygon(aes(fill = cell_type), color = "white") +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "right") +
    labs(fill = "Cell Type") +
    scale_fill_manual(values = cell_type_color_map)

  pdf(paste("./results/piechart/", names(df_list[i]),"_piechart_20percent.pdf",sep=""))
  print(pie_chart_plot)
  dev.off()

}
########################

##check heterogeneity
heterogeneity <- df_rescaled
heterogeneity$cluster <- NULL
heterogeneity$x_coord <- NULL
heterogeneity$y_coord <- NULL
heterogeneity$sample <- NULL

# Create a function to count non-zero values in a row
count_non_zero <- function(row) {
  return(sum(row > 0))
}

heterogeneity$number_cells <- apply(heterogeneity, 1, count_non_zero)

print(heterogeneity)

y <- se
y@meta.data[["heterogeneity"]] <- as.vector(heterogeneity$number_cells)

library(RColorBrewer)
color <- brewer.pal(4,"Spectral")
color <- rev(color)

pdf(file.path("./results/ST/gradient/heterogeneity_celltypes.pdf"))
print(FeatureOverlay(y, features = c("heterogeneity"), sampleids = 1:6,pt.size = 0.7,ncol = 2 , 
                   value.scale = "all" ,cols = color))
dev.off()



