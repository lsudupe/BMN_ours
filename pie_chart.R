


library(tidyverse)



se <- readRDS("./objects/heterogeneity/se_hierarchical.rds")

########cell types driving the cluster
x <- se

###coordinates
coor <- se@tools[["Staffli"]]@meta.data

##replace all values below 10% with 0
meta <- x@meta.data
a <- meta[,12:21]
a["cluster"] <- as.vector(x@meta.data[["clustering"]])

df <- a
# Filter values below 10%

# Calculate the 75th percentile threshold for each row
#thresholds <- apply(df[, -ncol(df)], 1, function(x) quantile(x, 0.75))

# Filter values below the 75th percentile threshold for each row
#df_filtered <- df %>%
 # mutate(across(-cluster, ~ifelse(. < thresholds[row_number()], 0, .)))

df["x_coord"] <- as.vector(coor$x) 
df["y_coord"] <- as.vector(coor$y) 
df["sample"] <- as.vector(meta$name) 

# Filter values below the 10% threshold for each row
df_filtered <- df %>%
  mutate(across(-c(cluster, x_coord, y_coord, sample), ~ifelse(. < 0.1, 0, .)))

# Rescale non-zero values in each row to sum up to 100% and overwrite the original values
#df_rescaled <- df_filtered %>%
 # mutate(across(-c(cluster, x_coord, y_coord, sample), ~ifelse(. != 0, . / sum(.[. != 0]) * 100, 0)))

# Rescale non-zero values in each row to sum up to 1 and overwrite the original values
df_rescaled <- df_filtered %>%
  rowwise() %>%
  mutate(across(-c(cluster, x_coord, y_coord, sample), 
                ~ifelse(. != 0, . / sum(c_across(-c(cluster, x_coord, y_coord, sample))[c_across(-c(cluster, x_coord, y_coord, sample)) != 0]), 0))) %>%
  ungroup()


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
# Reshape the data frame to a long format for plotting
data <- df_M8_F2_1C
data$sample <- NULL
data$cluster <- NULL
colnames(data)[which(names(data) == "x_coord")] <- "x"
colnames(data)[which(names(data) == "y_coord")] <- "y"

num_cell_types <- ncol(data) - 2
#data[, 1:num_cell_types] <- data[, 1:num_cell_types] / rowSums(data[, 1:num_cell_types])

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
radius <- 0.73  # Adjust this value to control the size of the pie charts

data_pie <- tidyr::expand_grid(x = unique(data_polar$x),
                               y = unique(data_polar$y),
                               point = seq_len(num_points + 1)) %>%
  left_join(data_polar, by = c("x", "y")) %>%
  mutate(angle = ifelse(point == 1, 0, start + (end - start) * (point - 2) / (num_points - 1)),
         x_plot = x + radius * cos(angle) * (point != 1),
         y_plot = y + radius * sin(angle) * (point != 1))

pie_chart_plot <- ggplot(data_pie, aes(x = x_plot, y = y_plot, group = interaction(x, y, cell_type))) +
  geom_polygon(aes(fill = cell_type), color = "white") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "right") +
  labs(fill = "Cell Type")

print(pie_chart_plot)
########################



##Ali approach
df_M1_fem_1C$sample <- NULL
plotting.data = cbind(df_M1_fem_1C, df_M1_fem_1C[12:13])
ggplot(plotting.data, aes(y=y_coord, x=x_coord, col=cluster))+
  geom_point()+
  theme_classic()

##Ali approach
plotting.data = cbind(B1@meta.data, B1@images$V11B18_363_B1@coordinates)
ggplot(plotting.data, aes(y=imagerow, x=imagecol, col=seurat_clusters))+
  geom_point()+
  theme_classic()


